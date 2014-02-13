#!/usr/bin/env python

import sys, os
import maf_parser as maf
import breakpoint_graph as bg
import permutations as perm

def get_synteny(maf_file, out_dir, min_block):
    out_permutations = os.path.join(out_dir, "genomes_permutations.txt")
    out_coords = os.path.join(out_dir, "blocks_coords.txt")
    #out_alignemnt = os.path.join(out_dir, "alignemnt_blocks.txt")

    MIN_ALIGNMENT = 100
    MAX_GAP = 5000

    permutations = maf.maf_to_permutations(maf_file, MIN_ALIGNMENT)
    lcb = bg.get_lcb(permutations, MAX_GAP)
    perm.filter_by_size(lcb, min_block)

    #perm.output_permutations(permutations, open(out_alignemnt, "w"))
    perm.output_permutations(lcb, open(out_permutations, "w"))
    perm.output_blocks_coords(lcb, open(out_coords, "w"))


def main():
    if len(sys.argv) != 4:
        sys.stderr.write("Usage: maf2synteny.py maf_file out_dir min_block\n")
        return

    maf_file = sys.argv[1]
    out_dir = sys.argv[2]
    min_block = sys.argv[3]
    get_synteny(maf_file, out_dir, int(min_block))


if __name__ == "__main__":
    main()
