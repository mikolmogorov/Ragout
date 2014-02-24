#!/usr/bin/env python

from __future__ import print_function
import sys, os
import maf_parser as maf
import breakpoint_graph as bg
import graph_compress
import bulges_removal
import permutations as perm


def get_synteny(maf_file, out_dir, min_block_size):
    out_permutations = os.path.join(out_dir, "genomes_permutations.txt")
    out_coords = os.path.join(out_dir, "blocks_coords.txt")

    MIN_ALIGNMENT = 100
    MAX_GAP = 5000
    MIN_FLANK = 500

    permutations = maf.maf_to_permutations(maf_file, MIN_ALIGNMENT)
    #permutations = perm.load_permutations(maf_file)

    blocks, block_groups = process_graph(permutations, MAX_GAP)
    perm.filter_by_size(blocks, min_block_size, MIN_FLANK, block_groups)

    perm.output_permutations(blocks, open(out_permutations, "w"))
    perm.output_blocks_coords(blocks, open(out_coords, "w"))


def process_graph(permutations, max_gap):
    graph = bg.build_graph(permutations)

    update = True
    while update:
        paths = graph_compress.compress_graph(graph, max_gap)
        bulges = bulges_removal.remove_bulges(graph, max_gap)
        update = bool(paths + bulges)
        print("Iteration:", paths, "paths", bulges, "bulges", file=sys.stderr)
    return graph.get_permutations()


def main():
    if len(sys.argv) != 4:
        print("Usage: maf2synteny.py maf_file out_dir min_block", file=sys.stderr)
        return

    maf_file = sys.argv[1]
    out_dir = sys.argv[2]
    min_block = sys.argv[3]
    get_synteny(maf_file, out_dir, int(min_block))


if __name__ == "__main__":
    main()
