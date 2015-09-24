#!/usr/bin/env python2.7

#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
A diagnostic script for chromosome fusions
"""

from __future__ import print_function
import sys
from collections import defaultdict

REF_NAME = "C57B6J"

def do_job(links_file, target_perms, all_blocks):
    block_chrs = {}
    with open(all_blocks, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                tokens = line[1:].split(".", 1)
                genome_name, chr_name = tokens
            elif genome_name != REF_NAME:
                continue
            else:
                blocks_ids = map(int, line.split(" ")[:-1])
                for b in blocks_ids:
                    block_chrs[abs(b)] = chr_name

    target_chrs = defaultdict(list)
    with open(target_perms, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                seq_name = line[1:]
            else:
                blocks_ids = map(int, line.split(" ")[:-1])
                for b in blocks_ids:
                    block_chr = block_chrs.get(abs(b), None)
                    if (not block_chr or not target_chrs[seq_name] or
                        target_chrs[seq_name][-1] != block_chr):
                        target_chrs[seq_name].append(str(b) + ":" + block_chr)

    with open(links_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("--") or line.startswith("sequence"):
                continue
            if line[0] not in "+-":
                print("\n" + line + "\n")
                continue

            sign, contig_name = line[0], line.split()[0][1:]
            blocks = (target_chrs[contig_name] if sign == "+"
                      else target_chrs[contig_name][::-1])
            print(sign, contig_name, blocks)


def main():
    if len(sys.argv) != 4:
        print("Usage: chromosome-report.py links target_perms all_perms")
        return 1
    do_job(sys.argv[1], sys.argv[2], sys.argv[3])
    return 0


if __name__ == "__main__":
    main()

