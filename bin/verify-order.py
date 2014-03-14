#!/usr/bin/env python

from __future__ import print_function
import sys
from collections import namedtuple, defaultdict
from itertools import product

from utils.nucmer_parser import *

Scaffold = namedtuple("Scaffold", ["name", "contigs"])
Contig = namedtuple("Contig", ["name", "sign", "gap"])


def parse_contigs_order(filename):
    scaffolds = []
    for line in open(filename, "r"):
        if line.startswith(">"):
            scaffolds.append(Scaffold(line.strip()[1:], []))
        else:
            name = line.strip("\n").replace("=", "_") #fix for nucmer
            without_sign = name[1:].strip()
            sign = 1 if name[0] == "+" else -1
            scaffolds[-1].contigs.append(Contig(without_sign, sign, 0))
    return scaffolds


def gap_count(lst_1, lst_2):
    if not lst_1 or not lst_2:
        return 0

    gaps = sys.maxint
    for i, j in product(lst_1, lst_2):
        gaps = min(gaps, abs(i.index - j.index) - 1)
    return gaps


def agreement(increasing, lst_1, lst_2, chr_len_dict):
    if not lst_1 or not lst_2:
        return True

    for i, j in product(lst_1, lst_2):
        same_chr = i.chr == j.chr
        if not same_chr:
            continue

        chr_len = chr_len_dict[i.chr]
        lower = chr_len * 0.2
        higher = chr_len * 0.8
        over_zero = ((increasing and i.coord > higher and j.coord < lower) or
                    (not increasing and i.coord < lower and j.coord > higher))

        if ((j.index > i.index) == increasing and
            abs(i.coord - j.coord) < chr_len / 3) or over_zero:
            return True
    return False


def do_job(nucmer_coords, scaffolds_ord):
    alignment = parse_nucmer_coords(nucmer_coords)
    #alignment = join_collinear_alignments(alignment)
    alignment = filter_by_coverage(alignment)
    entry_ord, chr_len, contig_len = get_order(alignment)
    scaffolds = parse_contigs_order(scaffolds_ord)

    #TODO: check signs
    total_breaks = 0
    total_gaps = 0
    total_contigs = 0
    for s in scaffolds:
        print("\n>" + s.name)

        prev = None
        increasing = None
        breaks = []
        for contig in s.contigs:
            if prev:
                if increasing is not None:
                    if not agreement(increasing, prev, entry_ord[contig.name], chr_len):
                        increasing = None
                        breaks.append(contig.name)
                        total_breaks += 1
                        sys.stdout.write("$")
                        #print map(lambda p: p.index, prev)
                else:
                    if len(entry_ord[contig.name]) == 1 and len(prev) == 1:
                        increasing = entry_ord[contig.name][0].index > prev[0].index

            print("{0}\t{1}\t{2}".format(contig.name, contig_len[contig.name],
                                         map(str, entry_ord[contig.name])))
            total_contigs += 1

            if gap_count(prev, entry_ord[contig.name]) > 0:
                total_gaps += 1
            #total_gaps += gap_count(prev, entry_ord[contig.name])

            #only if this contig has alignments
            if entry_ord[contig.name]:
                prev = entry_ord[contig.name]

        print("\tmiss-ordered: ", len(breaks))
    print("\nTotal miss-ordered:", total_breaks)
    print("Total gaps:", total_gaps)
    print("Total contigs:", total_contigs)
    print("Total scaffolds:", len(scaffolds))


def main():
    if len(sys.argv) < 3:
        print("Usage: verify-order.py <nucmer_coords> <scaffolds_ord>")
        return
    do_job(sys.argv[1], sys.argv[2])


if __name__ == "__main__":
    main()
