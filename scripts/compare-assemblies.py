#!/usr/bin/env python2.7

#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This script counts number of rearrangements between two assemblies
thus evaluating 'agreement level' between them
"""

from __future__ import print_function
from collections import defaultdict
import sys
import os
import argparse

import networkx as nx

from utils.lastz_parser import parse_lastz_maf, run_lastz
from utils.common import filter_intersecting, filter_by_length


MIN_ALIGNMENT = 5000


def get_alignment(reference, target, overwrite):
    out_file = (os.path.basename(reference) + "_" +
                os.path.basename(target) + ".maf")

    if os.path.isfile(out_file) and not overwrite:
        print("Alignment file already exists, lastz run skipped")
    else:
        run_lastz(reference, target, out_file)
    alignment = parse_lastz_maf(out_file)

    return alignment


def get_blocks(reference, target, overwrite):
    alignment = get_alignment(reference, target, overwrite)
    alignment = filter_by_length(alignment, MIN_ALIGNMENT)
    alignment = filter_intersecting(alignment)
    aln_with_id = list(enumerate(alignment))
    #alignment = join_collinear(alignment)

    #enumerating reference blocks
    ref_seqs = defaultdict(list)
    for aln_id, aln in aln_with_id:
        ref_seqs[aln.ref_id].append((aln_id, aln))
    for seq_id in ref_seqs:
        ref_seqs[seq_id].sort(key=lambda arec: arec[1].ref_start)
        to_block = (lambda (a_id, aln): (int(a_id) + 1) *
                    (1 if aln.ref_start < aln.ref_end else -1))
        ref_seqs[seq_id] = list(map(to_block, ref_seqs[seq_id]))

    #enumerating query blocks
    qry_seqs = defaultdict(list)
    for aln_id, aln in aln_with_id:
        qry_seqs[aln.qry_id].append((aln_id, aln))
    for seq_id in qry_seqs:
        qry_seqs[seq_id].sort(key=lambda arec: arec[1].qry_start)
        to_block = (lambda (a_id, aln): (int(a_id) + 1) *
                    (1 if aln.qry_start < aln.qry_end else -1))
        qry_seqs[seq_id] = list(map(to_block, qry_seqs[seq_id]))

    return ref_seqs, qry_seqs


def output_blocks(blocks):
    for seq, bl in blocks.items():
        print(">{0}".format(seq))
        print(" ".join(map(str,bl)))


def count_discord_adj(ref_blocks, qry_blocks):
    #building breakpoint graph
    graph = nx.MultiGraph()
    for seq, blocks in ref_blocks.items():
        for block_1, block_2 in zip(blocks[:-1], blocks[1:]):
            graph.add_edge(-block_1, block_2, name=seq, color="blue")
    for seq, blocks in qry_blocks.items():
        for block_1, block_2 in zip(blocks[:-1], blocks[1:]):
            graph.add_edge(-block_1, block_2, name=seq, color="green")

    counter = 0
    for node in graph.nodes():
        if len(graph.neighbors(node)) > 1:
            counter += 1

    return counter


def main():
    parser = argparse.ArgumentParser(description="Compare two assemblies")
    parser.add_argument("assembly_1", metavar="assembly_1",
                        help="path to first assembly")
    parser.add_argument("assembly_2", metavar="assembly_2",
                        help="path to second assembly")
    parser.add_argument("--overwrite", action="store_const", metavar="overwrite",
                        dest="overwrite", default=False, const=True,
                        help="overwrite existing lastz alignment")
    args = parser.parse_args()

    ref_blocks, qry_blocks = get_blocks(args.assembly_1, args.assembly_2,
                                        args.overwrite)
    #output_blocks(ref_blocks)
    #output_blocks(qry_blocks)
    print(count_discord_adj(ref_blocks, qry_blocks))


if __name__ == "__main__":
    main()
