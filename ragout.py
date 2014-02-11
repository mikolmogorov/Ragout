#!/usr/bin/env python

##################################
#The main Ragout module
##################################

import os
import sys
import logging
import argparse

import source.overlap as ovlp
import source.scaffolder as scfldr
import source.sibelia_parser as sp
import source.merge_iters as merge
import source.breakpoint_graph as bg
import source.config_parser as cparser
import source.assembly_refine as asref
import source.installer as installer
import source.utils as utils
from source.phylogeny import Phylogeny
from source.debug import DebugConfig
from source.permutation import PermutationContainer


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
        sys.stderr.write("Error: output directory doesn`t exists\n")
        return

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


def install_deps():
    installer.install_deps(LIB_DIR)


def check_dependencies():
    if not utils.which(SIBELIA_EXEC):
        sys.stderr.write("Sibelia is not installed. Use option \"--install-deps\" to install it.\n")
        return False
    return True


def main():
    parser = argparse.ArgumentParser(description="A tool for assisted assembly using multiple references")

    params = parser.add_argument_group()
    params.add_argument("-c", action="store", metavar="config", dest="config",
                        help="configuration file")

    params.add_argument("-o", action="store", metavar="output_dir", dest="output_dir",
                        help="output directory")

    params.add_argument("-s", action="store_const", metavar="skip_sibelia",
                        dest="skip_sibelia", default=False, const=True,
                        help="skip Sibelia running step")

    params.add_argument("-g", action="store_const", metavar="assembly_refine",
                        dest="assembly_refine", default=False, const=True,
                        help="refine with the assembly graph")

    parser.add_argument("--install-deps", action="store_const", metavar="install_deps",
                        dest="install_deps", default=False, const=True,
                        help="install Ragout dependencies")

    args = parser.parse_args()

    if args.install_deps:
        install_deps()
        return

    if not check_dependencies():
        return

    if not args.config:
        parser.print_usage()
        sys.stderr.write("error: argument -c is required\n")
        return
    if not args.output_dir:
        parser.print_usage()
        sys.stderr.write("error: argument -o is required\n")
        return

    do_job(args.config, args.output_dir, args.skip_sibelia, args.assembly_refine)

if __name__ == "__main__":
    main()
