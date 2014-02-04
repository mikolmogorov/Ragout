#!/usr/bin/env python

import sys
import maf
import breakpoint_graph as bg

def main():
    MIN_BLOCK = 500
    #permutations = maf.maf_to_permutations(sys.argv[1], MIN_BLOCK)
    #maf.output_permutations(permutations, open(sys.argv[2], "w"))
    permutations = maf.load_permutations(sys.argv[1])
    lcb = bg.get_lcb(permutations)
    maf.output_permutations(lcb, open(sys.argv[2], "w"))

if __name__ == "__main__":
    main()
