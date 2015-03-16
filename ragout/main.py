#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
The main Ragout module. It defines top-level logic of the program
"""

import os
import sys
import shutil
import logging
import argparse

import ragout.assembly_graph.assembly_refine as asref
import ragout.assembly_graph.assembly_graph as asgraph
import ragout.breakpoint_graph.breakpoint_graph as bg
import ragout.scaffolder.scaffolder as scfldr
import ragout.scaffolder.merge_iters as merge
import ragout.scaffolder.output_generator as out_gen
import ragout.maf2synteny.maf2synteny as m2s
import ragout.overlap.overlap as overlap
import ragout.shared.config as config
from ragout.overlap.overlap import OverlapException
from ragout.phylogeny.phylogeny import Phylogeny, PhyloException
from ragout.breakpoint_graph.permutation import (PermutationContainer,
                                                 PermException)
from ragout.synteny_backend.synteny_backend import (SyntenyBackend,
                                                    BackendException)
from ragout.parsers.recipe_parser import parse_ragout_recipe, RecipeException
from ragout.parsers.fasta_parser import read_fasta_dict, FastaError
from ragout.shared.debug import DebugConfig

#register backends
import synteny_backend.sibelia
import synteny_backend.cactus
import synteny_backend.maf
import synteny_backend.hal

logger = logging.getLogger()
debugger = DebugConfig.get_instance()


def enable_logging(log_file, debug):
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


def check_extern_modules(backend):
    """
    Checks if all neccessary native modules are available
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


def run_ragout(args):
    """
    A wrapper to catch possible exceptions
    """
    try:
        run_unsafe(args)
    except (RecipeException, PhyloException, PermException,
            BackendException, OverlapException, FastaError) as e:
        logger.error("An error occured while running Ragout:")
        logger.error(e)
        return 1

    return 0


def run_unsafe(args):
    """
    Top-level logic of the program
    """
    out_log = os.path.join(args.out_dir, "ragout.log")
    out_links = os.path.join(args.out_dir, "scaffolds.links")
    out_scaffolds = os.path.join(args.out_dir, "scaffolds.fasta")
    out_overlap = os.path.join(args.out_dir, "contigs_overlap.dot")
    out_refined_links = os.path.join(args.out_dir, "scaffolds_refined.links")
    out_refined_scaffolds = os.path.join(args.out_dir,
                                         "scaffolds_refined.fasta")
    debug_root = os.path.join(args.out_dir, "debug")

    if not os.path.isdir(args.out_dir):
        os.mkdir(args.out_dir)
    if args.debug:
        debugger.set_debug_dir(debug_root)
        debugger.clear_debug_dir()

    enable_logging(out_log, args.debug)
    logger.info("Cooking Ragout...")
    check_extern_modules(args.synteny_backend)
    recipe = parse_ragout_recipe(args.recipe)

    #Setting synteny block sizes
    synteny_blocks = config.vals["blocks"][recipe["blocks"]]

    #Running backend to get synteny blocks
    all_backends = SyntenyBackend.get_available_backends()
    backend = all_backends[args.synteny_backend]
    perm_files = backend.make_permutations(recipe, synteny_blocks, args.out_dir,
                                           args.overwrite, args.threads)

    #phylogeny-related
    if "tree" in recipe:
        logger.info("Phylogeny is taken from the recipe")
        phylogeny = Phylogeny.from_newick(recipe["tree"])
    else:
        logger.info("Inferring phylogeny from synteny blocks data")
        small_blocks = min(synteny_blocks)
        perm_container = PermutationContainer(perm_files[small_blocks],
                                              recipe, False, False, None)
        phylogeny = Phylogeny.from_permutations(perm_container)
        logger.info(phylogeny.tree_string)

    #main loop
    last_scaffolds = None
    for block_size in synteny_blocks:
        logger.info("Running Ragout with the block size {0}".format(block_size))

        if args.debug:
            debug_dir = os.path.join(debug_root, str(block_size))
            debugger.set_debug_dir(debug_dir)

        conservative = last_scaffolds is None
        perm_container = PermutationContainer(perm_files[block_size],
                                              recipe, args.resolve_repeats,
                                              conservative, phylogeny)

        graph = bg.BreakpointGraph()
        graph.build_from(perm_container, recipe)

        adjacencies = graph.find_adjacencies(phylogeny)
        scaffolds = scfldr.get_scaffolds(adjacencies, perm_container)

        if last_scaffolds is not None:
            last_scaffolds = merge.merge(last_scaffolds, scaffolds)
        else:
            last_scaffolds = scaffolds

    if args.debug:
        debugger.set_debug_dir(debug_root)

    logger.info("Reading contigs file")
    target_fasta_file = backend.get_target_fasta()
    target_fasta_dict = read_fasta_dict(target_fasta_file)

    out_gen.output_links(last_scaffolds, out_links)
    out_gen.output_fasta(target_fasta_dict, last_scaffolds, out_scaffolds)

    if not args.no_refine:
        overlap.make_overlap_graph(target_fasta_file, out_overlap)
        refined_scaffolds = asref.refine_scaffolds(out_overlap, last_scaffolds,
                                                   target_fasta_dict)
        out_gen.output_links(refined_scaffolds, out_refined_links)
        out_gen.output_fasta(target_fasta_dict, refined_scaffolds,
                             out_refined_scaffolds)
        if args.debug:
            shutil.copy(out_overlap, debugger.debug_dir)
            out_colored_overlap = os.path.join(debugger.debug_dir,
                                               "colored_overlap.dot")
            asgraph.save_colored_insert_overlap_graph(out_overlap, last_scaffolds,
                                                      refined_scaffolds,
                                                      out_colored_overlap)
        os.remove(out_overlap)

    logger.info("Your Ragout is ready!")


def main():
    parser = argparse.ArgumentParser(description="A tool for assisted assembly"
                                                 " using multiple references",
                                     formatter_class= \
                                        argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("recipe", metavar="recipe_file",
                        help="path to recipe file")
    parser.add_argument("-o", "--outdir", dest="out_dir",
                        metavar="output_dir",
                        help="path to the working directory",
                        default="ragout-out")
    parser.add_argument("-s", "--synteny", dest="synteny_backend",
                        default="sibelia", choices=["sibelia", "cactus", "maf", "hal"],
                        help="backend for synteny block decomposition")
    parser.add_argument("--no-refine", action="store_true",
                        dest="no_refine", default=False,
                        help="disable refinement with assembly graph")
    parser.add_argument("--overwrite", action="store_true", default=False,
                        dest="overwrite",
                        help="overwrite existing synteny blocks")
    parser.add_argument("--repeats", action="store_true", default=False,
                        dest="resolve_repeats",
                        help="try to resolve repeats before constructing BG")
    parser.add_argument("--debug", action="store_true",
                        dest="debug", default=False,
                        help="enable debug output")
    parser.add_argument("-t", "--threads", dest="threads", type=int,
                        default=1, help="number of threads for synteny backend")
    parser.add_argument("--version", action="version", version="Ragout v1.0")
    args = parser.parse_args()

    return run_ragout(args)
