#!/usr/bin/env python

##################################
#The main Ragout module
##################################

import os
import sys
import logging
import argparse

import src.overlap as ovlp
import src.scaffolder as scfldr
import src.merge_iters as merge
import src.breakpoint_graph as bg
import src.config_parser as cparser
import src.assembly_refine as asref
from src.synteny.synteny_backend import SyntenyBackend
from src.phylogeny import Phylogeny
from src.debug import DebugConfig
from src.permutation import PermutationContainer

#register backends
import src.synteny.sibelia
import src.synteny.cactus

logger = logging.getLogger()


def enable_logging(log_file):
    log_formatter = logging.Formatter("[%(asctime)s] %(name)s: %(levelname)s: %(message)s", "%H:%M:%S")
    console_formatter = logging.Formatter("[%(asctime)s] %(levelname)s: %(message)s", "%H:%M:%S")

    console_log = logging.StreamHandler()
    console_log.setLevel(logging.INFO)
    console_log.setFormatter(console_formatter)

    file_handler = logging.FileHandler(log_file, mode="w")
    file_handler.setFormatter(log_formatter)

    logger.setLevel(logging.DEBUG)
    logger.addHandler(console_log)
    logger.addHandler(file_handler)


#top-level logic of program
def do_job(config_file, out_dir, backend, assembly_refine):
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

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
    backends[backend].make_permutations(config, out_dir)

    last_scaffolds = None
    for block_size in config.blocks:
        logger.info("Running Ragout with the block size {0}...".format(block_size))
        block_dir = os.path.join(out_dir, str(block_size))
        block_config = os.path.join(block_dir, "blocks.cfg")
        block_order = os.path.join(block_dir, "scaffolds.ord")

        debug_dir = os.path.join(block_dir, "debug")
        DebugConfig.get_writer().set_debug_dir(debug_dir)


        perm_container = PermutationContainer(block_config)
        graph = bg.BreakpointGraph()
        graph.build_from(perm_container, True)

        connections = graph.find_adjacencies(phylogeny)
        scaffolds = scfldr.get_scaffolds(connections, perm_container)
        scfldr.output_order(scaffolds, block_order)

        if last_scaffolds:
            last_scaffolds = merge.merge(last_scaffolds, scaffolds)
        else:
            last_scaffolds = scaffolds

    scfldr.output_order(last_scaffolds, out_order)
    scfldr.output_scaffolds(config.targets, last_scaffolds, out_scaffolds)

    if assembly_refine:
        ovlp.make_overlap_graph(config.targets, out_overlap)
        refined_scaffolds = asref.refine_contigs(out_overlap, last_scaffolds)
        scfldr.output_order(refined_scaffolds, out_refined_order)
        scfldr.output_scaffolds(config.targets, refined_scaffolds, out_refined_scaffolds)

    logger.info("Your Ragout is ready!")


def main():
    parser = argparse.ArgumentParser(description="A tool for assisted assembly using multiple references",
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("config", metavar="config_file",
                        help="path to the configuration file")
    parser.add_argument("-o", "--outdir", dest="output_dir",
                        help="path to the working directory", default="ragout-out")
    parser.add_argument("-s", "--synteny", dest="synteny_backend", default="sibelia",
                        help="which tool to use for synteny block decomposition.",
                        choices=["sibelia", "cactus"])
    parser.add_argument("-r", "--refine", action="store_const", metavar="assembly_refine",
                        dest="assembly_refine", default=False, const=True,
                        help="refine with the assembly graph")
    parser.add_argument("-v", "--version", action="version", version="Ragout v0.1b")
    args = parser.parse_args()

    backends = SyntenyBackend.get_available_backends()
    if args.synteny_backend not in backends:
        sys.stderr.write(args.synteny_backend +
                         " is not installed. Use \"bin/install-deps.py\" to install it.\n")
        return

    do_job(args.config, args.output_dir, args.synteny_backend, args.assembly_refine)

if __name__ == "__main__":
    main()
