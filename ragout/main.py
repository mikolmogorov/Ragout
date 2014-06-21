#!/usr/bin/env python

"""
The main Ragout module
"""

import os
import sys
import logging
import argparse

import overlap.overlap as ovlp
import assembly_graph.assembly_refine as asref
import assembly_graph.assembly_graph as asgraph
import breakpoint_graph.breakpoint_graph as bg
from breakpoint_graph.phylogeny import Phylogeny, PhyloException
from breakpoint_graph.permutation import PermutationContainer, PermException
import scaffolder.scaffolder as scfldr
import scaffolder.merge_iters as merge
from synteny_backend.synteny_backend import SyntenyBackend
from parsers.recipe_parser import parse_ragout_recipe, RecipeException
from parsers.fasta_parser import read_fasta_dict, FastaError
from shared.debug import DebugConfig

#register backends
import synteny_backend.sibelia
import synteny_backend.cactus
import synteny_backend.maf

logger = logging.getLogger()
debugger = DebugConfig.get_instance()


def enable_logging(log_file, debug):
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


def do_job(recipe_file, out_dir, backend, assembly_refine,
           overwrite, debug):
    """
    top-level logic of program
    """
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    out_log = os.path.join(out_dir, "ragout.log")
    out_order = os.path.join(out_dir, "scaffolds.ord")
    out_scaffolds = os.path.join(out_dir, "scaffolds.fasta")
    out_overlap = os.path.join(out_dir, "contigs_overlap.dot")
    out_colored_overlap = os.path.join(out_dir, "colored_contigs_overlap.dot")
    out_compress_overlap = os.path.join(out_dir, "compress_contigs_overlap.dot")
    out_refined_order = os.path.join(out_dir, "scaffolds_refined.ord")
    out_refined_scaffolds = os.path.join(out_dir, "scaffolds_refined.fasta")

    if debug:
        core_debug = os.path.join(out_dir, "debug")
        debugger.set_debug_dir(core_debug)

    enable_logging(out_log, debug)
    logger.info("Cooking Ragout...")

    backends = SyntenyBackend.get_available_backends()

    try:
        recipe = parse_ragout_recipe(recipe_file)
    except RecipeException as e:
        logger.error("Error parsing recipe")
        logger.error(e)
        return 1

    try:
        phylogeny = Phylogeny(recipe)
    except PhyloException as e:
        logger.error(e)
        return 1

    perm_files = backends[backend].make_permutations(recipe, out_dir, overwrite)
    if not perm_files:
        logger.error("There were problems with synteny backend, exiting.")
        return 1


    last_scaffolds = None
    for block_size in recipe["blocks"]:
        logger.info("Running Ragout with the block size {0}".format(block_size))

        if debug:
            debug_dir = os.path.join(core_debug, str(block_size))
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

        if debug:
            ord_path = os.path.join(debug_dir, "scaffolds.ord")
            scfldr.output_order(scaffolds, ord_path)

        if last_scaffolds:
            last_scaffolds = merge.merge(last_scaffolds, scaffolds)
        else:
            last_scaffolds = scaffolds

    target_fasta_file = recipe["genomes"][recipe["target"]]["fasta"]
    try:
        target_fasta_dict = read_fasta_dict(target_fasta_file)
    except FastaError as e:
        logger.error(e)
        return 1

    scfldr.output_order(last_scaffolds, out_order)
    scfldr.output_fasta(target_fasta_dict, last_scaffolds, out_scaffolds)

    if assembly_refine:
        if not ovlp.make_overlap_graph(target_fasta_file, out_overlap):
            logger.error("Error in overlap graph reconstruction, exiting")
            return 1
        refined_scaffolds = asref.refine_scaffolds(out_overlap, last_scaffolds,
                                                   target_fasta_dict)
        #asgraph.save_colored_insert_overlap_graph(out_overlap, last_scaffolds,
        #                                          refined_scaffolds,
        #                                          out_colored_overlap)
        scfldr.output_order(refined_scaffolds, out_refined_order)
        scfldr.output_fasta(target_fasta_dict, refined_scaffolds,
                            out_refined_scaffolds)

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
    parser.add_argument("--refine", action="store_const", metavar="assembly_refine",
                        dest="assembly_refine", default=False, const=True,
                        help="enable refinement with assembly graph")
    parser.add_argument("--overwrite", action="store_const",
                        dest="overwrite", default=False, const=True,
                        help="overwrite existing synteny blocks")
    parser.add_argument("--debug", action="store_const",
                        dest="debug", default=False, const=True,
                        help="enable debug output")
    parser.add_argument("--version", action="version", version="Ragout v0.2b")
    args = parser.parse_args()

    backends = SyntenyBackend.get_available_backends()
    if args.synteny_backend not in backends:
        sys.stderr.write("\"{0}\" is not installed. You can use provided scripts "
                         "to install it.\n".format(args.synteny_backend))
        return 1

    return do_job(args.recipe, args.output_dir, args.synteny_backend,
                  args.assembly_refine, args.overwrite, args.debug)

if __name__ == "__main__":
    main()
