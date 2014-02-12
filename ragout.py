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
import src.sibelia_parser as sp
import src.merge_iters as merge
import src.breakpoint_graph as bg
import src.config_parser as cparser
import src.assembly_refine as asref
import src.utils as utils
from src.phylogeny import Phylogeny
from src.debug import DebugConfig
from src.permutation import PermutationContainer


LIB_DIR = "lib"
SIBELIA_DIR = os.path.join(LIB_DIR, "Sibelia")
SIBELIA_EXEC = "Sibelia"

running_dir = os.path.dirname(os.path.realpath(__file__))
os.environ["PATH"] += os.pathsep + os.path.join(running_dir, SIBELIA_DIR)

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
def do_job(config_file, out_dir, skip_sibelia, assembly_refine):
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    config = cparser.parse_ragout_config(config_file)
    phylogeny = Phylogeny(config.tree)

    out_log = os.path.join(out_dir, "log.txt")
    out_order = os.path.join(out_dir, "scaffolds.ord")
    out_scaffolds = os.path.join(out_dir, "scaffolds.fasta")
    out_overlap = os.path.join(out_dir, "contigs_overlap.dot")
    out_refined_order = os.path.join(out_dir, "scaffolds_refined.ord")
    out_refined_scaffolds = os.path.join(out_dir, "scaffolds_refined.fasta")

    enable_logging(out_log)
    last_scaffolds = None

    logger.info("Cooking Ragout...")
    for block_size in config.blocks:
        logger.info("Running with the block size {0}...".format(block_size))

        block_dir = os.path.join(out_dir, str(block_size))
        if not os.path.isdir(block_dir):
            os.mkdir(block_dir)
        block_config = os.path.join(block_dir, "blocks.cfg")
        block_order = os.path.join(block_dir, "scaffolds.ord")

        debug_dir = os.path.join(block_dir, "debug")
        DebugConfig.get_writer().set_debug_dir(debug_dir)

        if not skip_sibelia:
            sp.make_permutations(config.references, config.targets,
                                             block_size, block_dir)

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
        MIN_OVERLAP = 33
        MAX_PATH_LEN = 6
        ovlp.make_overlap_graph(config.targets, out_overlap, MIN_OVERLAP)
        refined_scaffolds = asref.refine_contigs(out_overlap, last_scaffolds, MAX_PATH_LEN)
        scfldr.output_order(refined_scaffolds, out_refined_order)
        scfldr.output_scaffolds(config.targets, refined_scaffolds, out_refined_scaffolds)

    logger.info("Your Ragout is ready!")


def check_dependencies():
    if not utils.which(SIBELIA_EXEC):
        return False
    return True


def main():
    parser = argparse.ArgumentParser(description="A tool for assisted assembly using multiple references")

    parser.add_argument("config", metavar="config_file",
                        help="path to the configuration file")

    parser.add_argument("-o", "--outdir", dest="output_dir",
                        help="path to the working directory [default = ragout-out]", default="ragout-out")

    #for debugging
    parser.add_argument("-s", "--skip-sibelia", action="store_true", dest="skip_sibelia",
                        help=argparse.SUPPRESS)

    parser.add_argument("-r", "--refine", action="store_const", metavar="assembly_refine",
                        dest="assembly_refine", default=False, const=True,
                        help="refine with the assembly graph")

    parser.add_argument("-v", "--version", action="version", version="Ragout v0.1b")

    args = parser.parse_args()

    if not check_dependencies():
        sys.stderr.write("Sibelia is not installed. Use \"bin/install-deps.py\" to install it.\n")
        return

    do_job(args.config, args.output_dir, args.skip_sibelia, args.assembly_refine)

if __name__ == "__main__":
    main()
