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
from breakpoint_graph.phylogeny import Phylogeny
from breakpoint_graph.permutation import PermutationContainer
import scaffolder.scaffolder as scfldr
import scaffolder.merge_iters as merge
from synteny_backend.synteny_backend import SyntenyBackend
import parsers.config_parser as cparser
from shared.debug import DebugConfig

#register backends
import synteny_backend.sibelia
import synteny_backend.cactus

logger = logging.getLogger()
debugger = DebugConfig.get_instance()


def enable_logging(log_file):
    log_formatter = logging.Formatter("[%(asctime)s] %(name)s: %(levelname)s: "
                                      "%(message)s", "%H:%M:%S")
    console_formatter = logging.Formatter("[%(asctime)s] %(levelname)s: "
                                          "%(message)s", "%H:%M:%S")
    console_log = logging.StreamHandler()
    console_log.setLevel(logging.INFO)
    console_log.setFormatter(console_formatter)

    file_handler = logging.FileHandler(log_file, mode="w")
    file_handler.setFormatter(log_formatter)

    logger.setLevel(logging.DEBUG)
    logger.addHandler(console_log)
    logger.addHandler(file_handler)


#top-level logic of program
def do_job(config_file, out_dir, backend, assembly_refine,
           circular_refs, overwrite, debug):
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    if debug:
        core_debug = os.path.join(out_dir, "debug")
        debugger.set_debug_dir(core_debug)

    config = cparser.parse_ragout_config(config_file)
    phylogeny = Phylogeny(config.tree)

    out_log = os.path.join(out_dir, "ragout.log")
    out_order = os.path.join(out_dir, "scaffolds.ord")
    out_scaffolds = os.path.join(out_dir, "scaffolds.fasta")
    out_overlap = os.path.join(out_dir, "contigs_overlap.dot")
    out_refined_order = os.path.join(out_dir, "scaffolds_refined.ord")
    out_refined_scaffolds = os.path.join(out_dir, "scaffolds_refined.fasta")

    enable_logging(out_log)

    logger.info("Cooking Ragout...")

    backends = SyntenyBackend.get_available_backends()
    if not backends[backend].make_permutations(config, out_dir, overwrite):
        logger.error("There were problems with synteny backend, exiting.")
        return

    last_scaffolds = None
    for block_size in config.blocks:
        logger.info("Running Ragout with the block size {0}".format(block_size))
        block_dir = os.path.join(out_dir, str(block_size))
        block_config = os.path.join(block_dir, "blocks.cfg")
        block_order = os.path.join(block_dir, "scaffolds.ord")

        if debug:
            debug_dir = os.path.join(core_debug, str(block_size))
            debugger.set_debug_dir(debug_dir)

        perm_container = PermutationContainer(block_config)
        graph = bg.BreakpointGraph()
        graph.build_from(perm_container, circular_refs)

        connections = graph.find_adjacencies(phylogeny)
        scaffolds = scfldr.get_scaffolds(connections, perm_container)
        scfldr.output_order(scaffolds, block_order)

        if last_scaffolds:
            last_scaffolds = merge.merge(last_scaffolds, scaffolds)
        else:
            last_scaffolds = scaffolds

    scfldr.output_order(last_scaffolds, out_order)
    scfldr.output_fasta(config.targets, last_scaffolds, out_scaffolds)

    if assembly_refine:
        if not ovlp.make_overlap_graph(config.targets, out_overlap):
            logger.error("Error in overlap graph reconstruction, exiting")
            return

        refined_scaffolds = asref.refine_scaffolds(out_overlap, last_scaffolds)
        scfldr.output_order(refined_scaffolds, out_refined_order)
        scfldr.output_fasta(config.targets, refined_scaffolds,
                                out_refined_scaffolds)

    logger.info("Your Ragout is ready!")


def main():
    parser = argparse.ArgumentParser(description="A tool for assisted assembly"
                                                 " using multiple references",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("config", metavar="config_file",
                        help="path to the configuration file")
    parser.add_argument("-o", "--outdir", dest="output_dir",
                        help="path to the working directory", default="ragout-out")
    parser.add_argument("-s", "--synteny", dest="synteny_backend", default="sibelia",
                        help="which tool to use for synteny block decomposition.",
                        choices=["sibelia", "cactus"])
    parser.add_argument("--refine", action="store_const", metavar="assembly_refine",
                        dest="assembly_refine", default=False, const=True,
                        help="refine with the assembly graph")
    parser.add_argument("--circular", action="store_const",
                        dest="circular_refs", default=False, const=True,
                        help="treat input references as circular")
    parser.add_argument("--overwrite", action="store_const",
                        dest="overwrite", default=False, const=True,
                        help="overwrite existing Sibelia/Cactus results")
    parser.add_argument("--debug", action="store_const",
                        dest="debug", default=False, const=True,
                        help="enable debug output")
    parser.add_argument("--version", action="version", version="Ragout v0.2b")
    args = parser.parse_args()

    backends = SyntenyBackend.get_available_backends()
    if args.synteny_backend not in backends:
        sys.stderr.write(args.synteny_backend + " is not installed."
                         "Use \"scripts/install-deps.py\" to install it.\n")
        return

    do_job(args.config, args.output_dir, args.synteny_backend,
           args.assembly_refine, args.circular_refs, args.overwrite, args.debug)

if __name__ == "__main__":
    main()
