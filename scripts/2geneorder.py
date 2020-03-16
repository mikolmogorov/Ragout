#!/usr/bin/python

#(c) 2013-2015 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This script converts synteny blocks file from Ragout internal format
to one that MLGO recognizes
"""

from __future__ import print_function
from __future__ import absolute_import
import sys
from collections import defaultdict

def main():
    if len(sys.argv) != 2:
        print("Usage: 2geneorder.py genome_permutations")
        return 1

    filename = sys.argv[1]
    blocks = defaultdict(list)
    with open(filename, "r") as f:
        while True:
            name_str = f.readline().strip()
            if not name_str:
                break

            blocks_str = f.readline().strip()
            blocks[name_str[1:].split(".")[0]].append(blocks_str)

    for genome, seqs in blocks.items():
        print(">{0}".format(genome))
        for seq in seqs:
            print(seq)

    return 0


if __name__ == "__main__":
    main()
