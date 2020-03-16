#!/usr/bin/env python2.7

#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This script checks correctness of Ragout results
if a 'true' reference is available
"""

from __future__ import print_function
from __future__ import absolute_import
import sys
import os
import argparse
from collections import namedtuple
from itertools import product

from utils.nucmer_parser import parse_nucmer_coords
from utils.common import (filter_by_coverage, join_collinear,
                          group_by_chr, get_order)
from six.moves import map

Scaffold = namedtuple("Scaffold", ["name", "contigs"])
Contig = namedtuple("Contig", ["name", "sign"])


def parse_links_file(filename):
    scaffolds = []

    def add_contig(string):
        name = string.replace("=", "_")      #fix for nucmer
        without_sign = name[1:].strip()
        if "[" in without_sign:
            without_sign = without_sign[:without_sign.index("[")]
        sign = 1 if name[0] == "+" else -1
        scaffolds[-1].contigs.append(Contig(without_sign, sign))

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("--") or line.startswith("sequence"):
                continue

            if line[0] not in ["+", "-"]:
                scaffolds.append(Scaffold(line, []))
            else:
                left_cont = line.split()[0]
                add_contig(left_cont)
        #add_contig(right_cont)

    return scaffolds


def parse_ord_file(filename):
    scaffolds = []

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                scaffolds.append(Scaffold(line[1:], []))
            else:
                sign = 1 if line[0] == "+" else -1
                scaffolds[-1].contigs.append(Contig(line[1:], sign))

    return scaffolds


def parse_contigs_order(filename):
    _filename, ext = os.path.splitext(filename)
    if ext[1:] == "ord":
        return parse_ord_file(filename)
    else:
        return parse_links_file(filename)


def gap_count(lst_1, lst_2):
    if not lst_1 or not lst_2:
        return 0

    gaps = sys.maxsize
    for i, j in product(lst_1, lst_2):
        gaps = min(gaps, abs(i.index - j.index) - 1)
    return gaps


def agreement_ord(increasing, lst_1, lst_2, chr_len_dict):
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
            abs(i.coord - j.coord) < chr_len / 1) or over_zero:
            return True
    return False


def agreement_strands(lst_1, lst_2, increasing):
    incr_sign = 1 if increasing else -1
    for i, j in product(lst_1, lst_2):
        if i == j and i == incr_sign:
            return True
    return False


def do_job(nucmer_coords, scaffolds_ord):
    alignment = parse_nucmer_coords(nucmer_coords)
    alignment = join_collinear(alignment)
    alignment = filter_by_coverage(alignment, 0.45)
    entry_ord, chr_len, contig_len = get_order(alignment)
    scaffolds = parse_contigs_order(scaffolds_ord)

    total_breaks = 0
    total_gaps = 0
    total_contigs = 0
    for s in scaffolds:
        print("\n>" + s.name)

        prev_aln = []
        prev_strand = None
        increasing = None
        breaks = []
        for contig in s.contigs:
            miss_ord = False
            miss_strand = False

            #checking order
            if prev_aln:
                if increasing is not None:
                    if not agreement_ord(increasing, prev_aln,
                                         entry_ord[contig.name], chr_len):
                        increasing = None
                        breaks.append(contig.name)
                        total_breaks += 1
                        miss_ord = True
                elif len(entry_ord[contig.name]) == 1 and len(prev_aln) == 1:
                    increasing = (entry_ord[contig.name][0].index >
                                  prev_aln[0].index)

            #checking strand
            cur_strand = [h.sign * contig.sign for h in entry_ord[contig.name]]
            if not miss_ord and prev_strand and cur_strand:
                if not agreement_strands(prev_strand, cur_strand, increasing):
                    breaks.append(contig.name)
                    total_breaks += 1
                    miss_strand = True
                    increasing = None

            if gap_count(prev_aln, entry_ord[contig.name]) > 0:
                total_gaps += 1

            #only if this contig has alignments
            if entry_ord[contig.name]:
                prev_aln = entry_ord[contig.name]
                prev_strand = cur_strand

            #output
            sign = "+" if contig.sign > 0 else "-"
            pos_list = list(map(str, entry_ord[contig.name]))
            pos_list_str = (str(pos_list) if len(pos_list) < 5 else
                            str(pos_list[:5]) + "...")
            print("{0}{1}\t{2}\t{3}".format(sign, contig.name,
                                            contig_len[contig.name], pos_list_str),
                                            end="")
            print("\t<<<order" if miss_ord else "", end="")
            print("\t<<<strand" if miss_strand else "", end="")
            print("")
            total_contigs += 1
            ###
        print("\tmiss-ordered: ", len(breaks))

    print("\nTotal miss-ordered:", total_breaks)
    print("Total gaps:", total_gaps)
    print("Total contigs:", total_contigs)
    print("Total scaffolds:", len(scaffolds))


def main():
    descr = ("A verification script for Ragout. It requires a contigs "
            "alignment on \"true\" reference in nucmer coords format. "
            "Given an alignment and a \"links\" file output by Ragout "
            "the script finds missassembled contigs and calculates some "
            "statistics.")
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument("nucmer_coords", metavar="nucmer_coords",
                        help="path to contigs alignment on 'true' reference")
    parser.add_argument("links_file", metavar="links_file",
                        help="path to 'links' file output by Ragout")
    args = parser.parse_args()

    do_job(args.nucmer_coords, args.links_file)


if __name__ == "__main__":
    main()
