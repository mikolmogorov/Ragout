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

import ragout.overlap.overlap as ovlp
import ragout.assembly_graph.assembly_refine as asref
import ragout.assembly_graph.assembly_graph as asgraph
import ragout.breakpoint_graph.breakpoint_graph as bg
import ragout.scaffolder.scaffolder as scfldr
import ragout.scaffolder.merge_iters as merge
import ragout.scaffolder.output_generator as out_gen
import ragout.overlap.overlap as overlap
import ragout.maf2synteny.maf2synteny as m2s
from ragout.breakpoint_graph.phylogeny import Phylogeny, PhyloException
from ragout.breakpoint_graph.permutation import (PermutationContainer,
                                                 PermException)
from ragout.synteny_backend.synteny_backend import SyntenyBackend
from ragout.parsers.recipe_parser import parse_ragout_recipe, RecipeException
from ragout.parsers.fasta_parser import read_fasta_dict, FastaError
from ragout.shared.debug import DebugConfig

#register backends
import synteny_backend.sibelia
import synteny_backend.cactus
import synteny_backend.maf

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
        logger.error("\"{0}\" is not installed. You can use provided scripts "
                     "to install it".format(backend))
        return False

    if not m2s.check_binary():
        return False

    if not overlap.check_binary():
        return False

    return True


def do_job(recipe_file, out_dir, backend, assembly_refine,
           overwrite, debug):
    """
    Top-level logic of program
    """
    out_log = os.path.join(out_dir, "ragout.log")
    out_links = os.path.join(out_dir, "scaffolds.links")
    out_scaffolds = os.path.join(out_dir, "scaffolds.fasta")
    out_overlap = os.path.join(out_dir, "contigs_overlap.dot")
    out_refined_links = os.path.join(out_dir, "scaffolds_refined.links")
    out_refined_scaffolds = os.path.join(out_dir, "scaffolds_refined.fasta")
    debug_root = os.path.join(out_dir, "debug")

    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    if debug:
        debugger.set_debug_dir(debug_root)
        debugger.clear_debug_dir()

    enable_logging(out_log, debug)
    logger.info("Cooking Ragout...")

    if not check_extern_modules(backend):
        return 1

    try:
        recipe = parse_ragout_recipe(recipe_file)
        phylogeny = Phylogeny(recipe)
        target_fasta_file = recipe["genomes"][recipe["target"]]["fasta"]
        logger.info("Reading FASTA with contigs")
        target_fasta_dict = read_fasta_dict(target_fasta_file)

    except (RecipeException, FastaError, PhyloException) as e:
        logger.error(e)
        return 1

    backends = SyntenyBackend.get_available_backends()
    perm_files = backends[backend].make_permutations(recipe, out_dir, overwrite)
    if not perm_files:
        logger.error("There were problems with synteny backend, exiting.")
        return 1

    last_scaffolds = None
    for block_size in recipe["blocks"]:
        logger.info("Running Ragout with the block size {0}".format(block_size))

        if debug:
            debug_dir = os.path.join(debug_root, str(block_size))
            debugger.set_debug_dir(debug_dir)

        try:
            perm_container = PermutationContainer(perm_files[block_size],
                                                  recipe)
        except PermException as e:
            logger.error(e)
            return 1

        graph = bg.BreakpointGraph()
        graph.build_from(perm_container, recipe)

        adjacencies = graph.find_adjacencies(phylogeny)
        scaffolds = scfldr.get_scaffolds(adjacencies, perm_container)

        if last_scaffolds:
            last_scaffolds = merge.merge(last_scaffolds, scaffolds)
        else:
            last_scaffolds = scaffolds

    debugger.set_debug_dir(debug_root)
    out_gen.output_links(last_scaffolds, out_links)
    out_gen.output_fasta(target_fasta_dict, last_scaffolds, out_scaffolds)

    if assembly_refine:
        if not ovlp.make_overlap_graph(target_fasta_file, out_overlap):
            logger.error("Error in overlap graph reconstruction, exiting")
            return 1
        refined_scaffolds = asref.refine_scaffolds(out_overlap, last_scaffolds,
                                                   target_fasta_dict)
        out_gen.output_links(refined_scaffolds, out_refined_links)
        out_gen.output_fasta(target_fasta_dict, refined_scaffolds,
                            out_refined_scaffolds)
        if debug:
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
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("recipe", metavar="recipe_file",
                        help="path to recipe file")
    parser.add_argument("-o", "--outdir", dest="output_dir",
                        help="path to the working directory",
                        default="ragout-out")
    parser.add_argument("-s", "--synteny", dest="synteny_backend",
                        default="sibelia", choices=["sibelia", "cactus", "maf"],
                        help="which tool to use for synteny block decomposition")
    parser.add_argument("--no-refine", action="store_const", metavar="no_refine",
                        dest="no_refine", default=False, const=True,
                        help="disable refinement with assembly graph")
    parser.add_argument("--overwrite", action="store_const",
                        dest="overwrite", default=False, const=True,
                        help="overwrite existing synteny blocks")
    parser.add_argument("--debug", action="store_const",
                        dest="debug", default=False, const=True,
                        help="enable debug output")
    parser.add_argument("--version", action="version", version="Ragout v0.3b")
    args = parser.parse_args()

    return do_job(args.recipe, args.output_dir, args.synteny_backend,
                  not args.no_refine, args.overwrite, args.debug)
