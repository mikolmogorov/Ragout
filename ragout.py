#!/usr/bin/env python

import os
import argparse

import source.breakpoint_graph as bg
import source.sibelia_parser as sp
import source.overlap as ovlp
import source.debrujin_refine as debrujin
import source.scaffolder as scfldr
import source.merge_iters as merge
from source.phylogeny import Phylogeny
from source.permutation import PermutationContainer, parse_config
from source.debug import DebugConfig

SIBELIA_BIN = "/home/volrath/Bioinf/Sibelia/distr/bin/"
os.environ["PATH"] += os.pathsep + os.path.abspath(SIBELIA_BIN)

def do_job(config_file, out_dir, skip_sibelia, debrujin_refine):
    if not os.path.isdir(out_dir):
        sys.stderr.write("Output directory doesn`t exists\n")
        return

    references, targets, tree_string, block_sizes = parse_config(config_file)
    phylogeny = Phylogeny(tree_string)

    out_order = os.path.join(out_dir, "scaffolds.ord")
    out_scaffolds = os.path.join(out_dir, "scaffolds.fasta")
    out_overlap = os.path.join(out_dir, "contigs_overlap.dot")
    out_refined_order = os.path.join(out_dir, "scaffolds_refined.ord")
    out_refined_scaffolds = os.path.join(out_dir, "scaffolds_refined.fasta")
    oout_scaffolds = os.path.join(out_dir, "scaffolds.fasta")

    last_scaffolds = None

    for block_size in block_sizes:
        block_dir = os.path.join(out_dir, str(block_size))
        if not os.path.isdir(block_dir):
            os.mkdir(block_dir)
        block_config = os.path.join(block_dir, "blocks.cfg")
        block_order = os.path.join(block_dir, "scaffolds.ord")

        debug_dir = os.path.join(block_dir, "debug")
        DebugConfig.get_writer().set_debug_dir(debug_dir)

        if not skip_sibelia:
            sp.make_permutations(references, targets, block_size, block_dir)

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
    scfldr.output_scaffolds(targets, last_scaffolds, out_scaffolds)

    if debrujin_refine:
        ovlp.make_overlap_graph(targets, out_overlap)
        refined_scaffolds = debrujin.refine_contigs(out_overlap, last_scaffolds)
        scfldr.output_order(refined_scaffolds, out_refined_order)
        scfldr.output_scaffolds(targets, refined_scaffolds, out_refined_scaffolds)


def main():
    parser = argparse.ArgumentParser(description="Tool for reference-assisted assembly")
    parser.add_argument("-c", action="store", metavar="config", dest="config",
                        required=True, help="Configuration file")
    parser.add_argument("-o", action="store", metavar="output_dir", dest="output_dir",
                        required=True, help="Output directory")
    parser.add_argument("-s", action="store_const", metavar="skip_sibelia", dest="skip_sibelia",
                        default=False, const=True, help="Skip Sibelia running step")
    parser.add_argument("-g", action="store_const", metavar="debrujin_refine", dest="debrujin_refine",
                        default=False, const=True, help="Refine with Debrujin graph")

    args = parser.parse_args()
    do_job(args.config, args.output_dir, args.skip_sibelia, args.debrujin_refine)

if __name__ == "__main__":
    main()
