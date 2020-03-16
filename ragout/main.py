#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
The main Ragout module. It defines top-level logic of the program
"""

from __future__ import absolute_import
from __future__ import division
import os
import shutil
import logging
import argparse
from collections import namedtuple

import ragout.assembly_graph.assembly_refine as asref
import ragout.scaffolder.scaffolder as scfldr
import ragout.scaffolder.merge_iters as merge
import ragout.maf2synteny.maf2synteny as m2s
import ragout.overlap.overlap as overlap
import ragout.shared.config as config
from ragout.scaffolder.output_generator import OutputGenerator
from ragout.overlap.overlap import OverlapException
from ragout.phylogeny.phylogeny import Phylogeny
from ragout.parsers.phylogeny_parser import PhyloException
from ragout.breakpoint_graph.permutation import (PermutationContainer,
                                                 PermException)
from ragout.synteny_backend.synteny_backend import (SyntenyBackend,
                                                    BackendException)
from ragout.parsers.recipe_parser import parse_ragout_recipe, RecipeException
from ragout.parsers.fasta_parser import read_fasta_dict, FastaError
from ragout.shared.debug import DebugConfig
from ragout.breakpoint_graph.breakpoint_graph import BreakpointGraph
from ragout.breakpoint_graph.inferer import AdjacencyInferer
from ragout.breakpoint_graph.chimera_detector import ChimeraDetector
from ragout.__version__ import __version__

#register backends
import ragout.synteny_backend.sibelia
#import synteny_backend.cactus
import ragout.synteny_backend.maf
import ragout.synteny_backend.hal
import ragout.six as six

logger = logging.getLogger()
debugger = DebugConfig.get_instance()


RunStage = namedtuple("RunStage", ["name", "block_size", "ref_indels",
                                   "repeats", "rearrange"])
ID_SMALLEST = -1


def _enable_logging(log_file, debug):
    """
    Turns on logging, sets debug levels and assigns a log file
    """
    log_formatter = logging.Formatter("[%(asctime)s] %(name)s: %(levelname)s: "
                                      "%(message)s", "%H:%M:%S")
    console_formatter = logging.Formatter("[%(asctime)s] %(levelname)s: "
                                          "%(message)s", "%H:%M:%S")
    console_log = logging.StreamHandler()
    console_log.setFormatter(console_formatter)
    if not debug:
        console_log.setLevel(logging.INFO)

    file_handler = logging.FileHandler(log_file, mode="w")
    file_handler.setFormatter(log_formatter)

    logger.setLevel(logging.DEBUG)
    logger.addHandler(console_log)
    logger.addHandler(file_handler)


def _check_extern_modules(backend):
    """
    Checks if all necessary native modules are available
    """
    backends = SyntenyBackend.get_available_backends()
    if backend not in backends:
        raise BackendException("\"{0}\" is not installed.".format(backend))

    if not m2s.check_binary():
        raise BackendException("maf2synteny binary is missing, "
                               "did you run 'make'?")

    if not overlap.check_binary():
        raise BackendException("overlap binary is missing, "
                               "did you run 'make'?")


def _make_run_stages(block_sizes, resolve_repeats):
    """
    Setting parameters of run stages (iterations)
    """
    stages = []
    #general stages for structural assembly
    for block in block_sizes:
        stages.append(RunStage(name=str(block), block_size=block,
                               ref_indels=False, repeats=False,
                               rearrange=True))

    #refining stage to close assembly gaps
    stages.append(RunStage(name="refine", block_size=block_sizes[ID_SMALLEST],
                           ref_indels=True, repeats=resolve_repeats,
                           rearrange=False))
    return stages


def _get_phylogeny_and_naming_ref(recipe, permutation_file):
    """
    Retrieves phylogeny (infers if necessary) as well as
    naming reference genome
    """
    if "tree" in recipe:
        logger.info("Phylogeny is taken from the recipe")
        phylogeny = Phylogeny.from_newick(recipe["tree"])
    else:
        logger.info("Inferring phylogeny from synteny blocks data")
        perm_cont = PermutationContainer(permutation_file, recipe,
                                         resolve_repeats=False,
                                         allow_ref_indels=True,
                                         phylogeny=None)
        phylogeny = Phylogeny.from_permutations(perm_cont)
        logger.info("Inferred tree: %s", phylogeny.tree_string)

    if "naming_ref" in recipe:
        naming_ref = recipe["naming_ref"]
    else:
        leaves_sorted = phylogeny.leaves_by_distance(recipe["target"])
        naming_ref = leaves_sorted[0]
        logger.info("'%s' is chosen as a naming reference", naming_ref)

    return phylogeny, naming_ref


def _get_synteny_scale(recipe, synteny_backend):
    if "blocks" in recipe:
        if isinstance(recipe["blocks"], six.string_types):
            scale = config.vals["blocks"][recipe["blocks"]]
        else:
            scale = recipe["blocks"]
    else:
        scale = config.vals["blocks"][synteny_backend.infer_block_scale(recipe)]

    logger.info("Running withs synteny block sizes '%s'", str(scale))
    return scale


def _run_ragout(args):
    """
    Top-level logic of the program
    """
    if not os.path.isdir(args.out_dir):
        os.mkdir(args.out_dir)

    debug_root = os.path.join(args.out_dir, "debug")
    debugger.set_debugging(args.debug)
    debugger.set_debug_dir(debug_root)
    debugger.clear_debug_dir()

    out_log = os.path.join(args.out_dir, "ragout.log")
    _enable_logging(out_log, args.debug)
    logger.info("Starting Ragout v%s", str(__version__))

    #parsing recipe, preparing synteny backend
    _check_extern_modules(args.synteny_backend)
    synteny_backend = SyntenyBackend.get_available_backends() \
                                        [args.synteny_backend]
    recipe = parse_ragout_recipe(args.recipe)
    synteny_sizes = _get_synteny_scale(recipe, synteny_backend)

    #Running synteny backend to get synteny blocks
    perm_files = synteny_backend.make_permutations(recipe, synteny_sizes,
                                                   args.out_dir, args.overwrite,
                                                   args.threads)

    #setting up phylogenetic tree
    phylo_perm_file = perm_files[synteny_sizes[ID_SMALLEST]]
    phylogeny, naming_ref = _get_phylogeny_and_naming_ref(recipe,
                                                          phylo_perm_file)

    #parsing permutation files, apply filters and build breakpoint graph
    logger.info("Processing permutation files")
    raw_bp_graphs = {}
    stage_perms = {}
    run_stages = _make_run_stages(synteny_sizes, args.resolve_repeats)
    for stage in run_stages:
        debugger.set_debug_dir(os.path.join(debug_root, stage.name))
        stage_perms[stage] = PermutationContainer(perm_files[stage.block_size],
                                                  recipe, stage.repeats,
                                                  stage.ref_indels, phylogeny)
        raw_bp_graphs[stage] = BreakpointGraph(stage_perms[stage])

    #initializing chimera detector
    target_sequences = read_fasta_dict(synteny_backend.get_target_fasta())
    chim_detect = None
    if not args.solid_scaffolds:
        chim_detect = ChimeraDetector(raw_bp_graphs, run_stages, target_sequences)

    #inferring adjacencies: loop over stages (iterations)
    scaffolds = None
    prev_stages = []
    for stage in run_stages:
        logger.info("Stage \"%s\"", stage.name)
        debugger.set_debug_dir(os.path.join(debug_root, stage.name))
        prev_stages.append(stage)

        broken_perms = stage_perms[stage]
        if not args.solid_scaffolds:
            broken_perms = chim_detect.break_contigs(stage_perms[stage], [stage])
        breakpoint_graph = BreakpointGraph(broken_perms)

        adj_inferer = AdjacencyInferer(breakpoint_graph, phylogeny)
        adjacencies = adj_inferer.infer_adjacencies()
        cur_scaffolds = scfldr.build_scaffolds(adjacencies, broken_perms,
                                               debug_output=True,
                                               correct_distances=False)

        if scaffolds is not None:
            if not args.solid_scaffolds:
                broken_perms = chim_detect.break_contigs(stage_perms[stage],
                                                         prev_stages)
            cur_scaffolds = merge.merge_scaffolds(scaffolds, cur_scaffolds,
                                                  broken_perms,
                                                  stage.rearrange)

        merge.get_breakpoints(cur_scaffolds, breakpoint_graph, broken_perms)

        scaffolds = cur_scaffolds

    debugger.set_debug_dir(debug_root)
    ####

    #name scaffolds according to one of the references
    last_stage = run_stages[ID_SMALLEST]
    scfldr.assign_scaffold_names(scaffolds, stage_perms[last_stage], naming_ref)
    scfldr.update_gaps(scaffolds)

    #refine with the assembly graph
    if args.refine:
        out_overlap = os.path.join(args.out_dir, "contigs_overlap.dot")
        overlap.make_overlap_graph(synteny_backend.get_target_fasta(),
                                   out_overlap)
        scaffolds = asref.refine_scaffolds(out_overlap, scaffolds,
                                           target_sequences)
        if args.debug:
            shutil.copy(out_overlap, debugger.debug_dir)
        os.remove(out_overlap)

    out_gen = OutputGenerator(target_sequences, scaffolds)
    out_gen.make_output(args.out_dir, recipe["target"])
    logger.info("Done!")


def main():
    parser = argparse.ArgumentParser(
                description="Chromosome assembly with multiple "
                "references",
                formatter_class= argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("recipe", metavar="recipe_file",
                        help="path to recipe file")
    parser.add_argument("-o", "--outdir", dest="out_dir",
                        metavar="output_dir",
                        help="output directory",
                        default="ragout-out")
    parser.add_argument("-s", "--synteny", dest="synteny_backend",
                        default="sibelia",
                        choices=["sibelia", "maf", "hal"],
                        help="backend for synteny block decomposition")
    parser.add_argument("--refine", action="store_true",
                        dest="refine", default=False,
                        help="enable refinement with assembly graph")
    parser.add_argument("--solid-scaffolds", action="store_true",
                        dest="solid_scaffolds", default=False,
                        help="do not break input sequences - disables "
                        "chimera detection module")
    parser.add_argument("--overwrite", action="store_true", default=False,
                        dest="overwrite",
                        help="overwrite results from the previous run")
    parser.add_argument("--repeats", action="store_true", default=False,
                        dest="resolve_repeats",
                        help="enable repeat resolution algorithm")
    parser.add_argument("--debug", action="store_true",
                        dest="debug", default=False,
                        help="enable debug output")
    parser.add_argument("-t", "--threads", dest="threads", type=int,
                        default=1, help="number of threads for synteny backend")
    parser.add_argument("--version", action="version", version=__version__)
    args = parser.parse_args()

    try:
        _run_ragout(args)
    except (RecipeException, PhyloException, PermException,
            BackendException, OverlapException, FastaError) as e:
        logger.error("An error occured while running Ragout:")
        logger.error(e)
        return 1

    return 0
