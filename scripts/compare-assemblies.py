#!/usr/bin/env python2.7

#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This script counts number of rearrangements between two assemblies
thus evaluating 'agreement level' between them
"""

from __future__ import print_function
from __future__ import absolute_import
from collections import defaultdict
import sys
import os
import argparse

import networkx as nx

from utils.lastz_parser import (parse_lastz_maf, run_lastz,
                                filter_intersecting, filter_by_length)
from six.moves import map
from six.moves import zip


def get_alignment(reference, target, overwrite):
    out_file = (os.path.basename(reference) + "_" +
                os.path.basename(target) + ".maf")

    if os.path.isfile(out_file) and not overwrite:
        print("Alignment file already exists, lastz run skipped")
    else:
        run_lastz(reference, target, out_file)
    alignment = parse_lastz_maf(out_file)

    return alignment


def get_blocks(reference, target, overwrite, min_alignmtnt):
    alignment = get_alignment(reference, target, overwrite)
    alignment = filter_by_length(alignment, min_alignmtnt)
    alignment = filter_intersecting(alignment)
    #alignment = join_collinear(alignment)

    def enum_blocks(aln_rows):
        blocks = defaultdict(list)
        for r_id, row in enumerate(aln_rows):
            blocks[row.seq_id].append((r_id, row))
        for seq_id in blocks:
            blocks[seq_id].sort(key=lambda pair: pair[1].start)
            to_block = lambda r_id_row: (r_id_row[0] + 1) * r_id_row[1].strand
            blocks[seq_id] = list(map(to_block, blocks[seq_id]))

        return blocks

    #IMPORTANT: ref/qry rows should have corresponding order
    ref_blocks = enum_blocks([ap.ref for ap in alignment])
    qry_blocks = enum_blocks([ap.qry for ap in alignment])

    return ref_blocks, qry_blocks


def output_blocks(blocks):
    for seq, bl in blocks.items():
        print(">{0}".format(seq))
        print(" ".join(["{0:+d}".format(b) for b in bl]))


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
    for node in graph.nodes:
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
    parser.add_argument("-b", "--block", dest="block_size",
                        help="minimum synteny block size",
                        default="5000")
    args = parser.parse_args()

    ref_blocks, qry_blocks = get_blocks(args.assembly_1, args.assembly_2,
                                        args.overwrite, int(args.block_size))
    #output_blocks(ref_blocks)
    #output_blocks(qry_blocks)
    print(count_discord_adj(ref_blocks, qry_blocks))


if __name__ == "__main__":
    main()
