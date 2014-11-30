#!/usr/bin/env python2.7

#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This script counts number of rearrangements between two assemblies
thus evaluating 'agreement level' between them
"""

from __future__ import print_function
import sys

from utils.nucmer_parser import *
import networkx as nx


def get_blocks(nucmer_coords):
    def sign(coord_1, coord_2):
        return 1 if coord_1 < coord_2 else -1

    #TODO: filter repeats
    #TODO: chech they don't intersect
    alignment = parse_nucmer_coords(nucmer_coords)
    #alignment = join_collinear(alignment)
    aln_with_id = list(enumerate(filter_by_coverage(alignment)))

    #enumerating reference blocks
    ref_seqs = defaultdict(list)
    for aln_id, aln in aln_with_id:
        ref_seqs[aln.ref_id].append((aln_id, aln))
    for seq_id in ref_seqs:
        ref_seqs[seq_id].sort(key=lambda arec: arec[1].s_ref)
        to_block = (lambda (a_id, aln): (int(a_id) + 1) *
                    (1 if aln.s_ref < aln.e_ref else -1))
        ref_seqs[seq_id] = list(map(to_block, ref_seqs[seq_id]))

    #enumerating query blocks
    qry_seqs = defaultdict(list)
    for aln_id, aln in aln_with_id:
        qry_seqs[aln.qry_id].append((aln_id, aln))
    for seq_id in qry_seqs:
        qry_seqs[seq_id].sort(key=lambda arec: arec[1].s_qry)
        to_block = (lambda (a_id, aln): (int(a_id) + 1) *
                    (1 if aln.s_qry < aln.e_qry else -1))
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
            #print(node, graph.neighbors(node))
            #for n in graph.neighbors(node):
            #    print(graph[node][n])

    return counter


def main():
    ref_blocks, qry_blocks = get_blocks(sys.argv[1])
    #output_blocks(ref_blocks)
    #output_blocks(qry_blocks)
    print(count_discord_adj(ref_blocks, qry_blocks))


if __name__ == "__main__":
    main()
