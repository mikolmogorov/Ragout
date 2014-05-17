#!/usr/bin/env python

##################################
#The main Ragout module
##################################

import os
import sys
import logging
import argparse

import overlap.overlap as ovlp
import assembly_graph.assembly_refine as asref
import breakpoint_graph.breakpoint_graph as bg
from breakpoint_graph.phylogeny import Phylogeny, PhyloException
from breakpoint_graph.permutation import PermutationContainer, PermException
import scaffolder.scaffolder as scfldr
import scaffolder.merge_iters as merge
from synteny_backend.synteny_backend import SyntenyBackend
from parsers.recipe_parser import parse_ragout_recipe, RecipeException
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


#top-level logic of program
def do_job(recipe_file, out_dir, backend, assembly_refine,
           circular_refs, overwrite, debug):
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    out_log = os.path.join(out_dir, "ragout.log")
    out_order = os.path.join(out_dir, "scaffolds.ord")
    out_scaffolds = os.path.join(out_dir, "scaffolds.fasta")
    out_overlap = os.path.join(out_dir, "contigs_overlap.dot")
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
        return
    #if recipe.maf:
    #    backend = "maf"

    try:
        phylogeny = Phylogeny(recipe)
    except PhyloException as e:
        logger.error(e)
        return

    perm_files = backends[backend].make_permutations(recipe, out_dir, overwrite)
    if not perm_files:
        logger.error("There were problems with synteny backend, exiting.")
        return

    last_scaffolds = None
    for block_size in recipe.blocks:
        logger.info("Running Ragout with the block size {0}".format(block_size))

        if debug:
            debug_dir = os.path.join(core_debug, str(block_size))
            debugger.set_debug_dir(debug_dir)

        try:
            perm_container = PermutationContainer(perm_files[block_size],
                                                  recipe)
        except PermException as e:
            logger.error(e)
            return

        graph = bg.BreakpointGraph()
        graph.build_from(perm_container, circular_refs)

        connections = graph.find_adjacencies(phylogeny)
        scaffolds = scfldr.get_scaffolds(connections, perm_container)

        if debug:
            ord_path = os.path.join(debug_dir, "scaffolds.ord")
            scfldr.output_order(scaffolds, ord_path)

        if last_scaffolds:
            last_scaffolds = merge.merge(last_scaffolds, scaffolds)
        else:
            last_scaffolds = scaffolds

    target_fasta = recipe.fasta[recipe.targets[0]]
    scfldr.output_order(last_scaffolds, out_order)
    scfldr.output_fasta(target_fasta, last_scaffolds, out_scaffolds)

    if assembly_refine:
        if not ovlp.make_overlap_graph(target_fasta, out_overlap):
            logger.error("Error in overlap graph reconstruction, exiting")
            return

        refined_scaffolds = asref.refine_scaffolds(out_overlap, last_scaffolds)
        scfldr.output_order(refined_scaffolds, out_refined_order)
        scfldr.output_fasta(target_fasta, refined_scaffolds,
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
                        help="tool for synteny block decomposition")
    parser.add_argument("--refine", action="store_const", metavar="assembly_refine",
                        dest="assembly_refine", default=False, const=True,
                        help="enable refinement with assembly graph")
    parser.add_argument("--circular", action="store_const",
                        dest="circular_refs", default=False, const=True,
                        help="treat input references as circular")
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
        sys.stderr.write(args.synteny_backend + " is not installed. "
                         "You can use provided scripts to install it.\n")
        return

    do_job(args.recipe, args.output_dir, args.synteny_backend,
           args.assembly_refine, args.circular_refs, args.overwrite, args.debug)

if __name__ == "__main__":
    main()
