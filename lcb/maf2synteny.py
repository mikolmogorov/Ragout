#!/usr/bin/env python

import sys, os
import maf_parser as maf
import breakpoint_graph as bg
import permutations as perm

def main():
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: maf2synteny.py maf_file out_dir\n")
        return

    maf_file = sys.argv[1]
    out_dir = sys.argv[2]
    out_permutations = os.path.join(out_dir, "genomes_permutations.txt")
    out_coords = os.path.join(out_dir, "blocks_coords.txt")
    out_alignemnt = os.path.join(out_dir, "alignemnt_blocks.txt")

    MIN_ALIGNMENT = 100
    MIN_BLOCK = 500

    permutations = maf.maf_to_permutations(maf_file, MIN_ALIGNMENT)
    #perm.output_permutations(permutations, open(out_alignemnt, "w"))

    lcb = bg.get_lcb(permutations, MIN_BLOCK)

    perm.output_permutations(lcb, open(out_permutations, "w"))
    perm.output_blocks_coords(lcb, open(out_coords, "w"))

if __name__ == "__main__":
    main()
