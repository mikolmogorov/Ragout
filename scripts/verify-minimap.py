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
from collections import namedtuple
from itertools import product
import subprocess

from utils.common import AlignmentRow, AlignmentColumn
from utils.common import (filter_by_coverage, join_collinear,
                          group_by_chr, get_order)
from six.moves import map
from six.moves import range

Scaffold = namedtuple("Scaffold", ["name", "contigs"])
Contig = namedtuple("Contig", ["name", "sign"])


MINIMAP_BIN = "minimap2"


def read_fasta_dict(filename):
    """
    Reads fasta file into dictionary. Also preforms some validation
    """
    #logger.debug("Reading contigs file")

    header = None
    seq = []
    fasta_dict = {}

    try:
        with open(filename, "r") as f:
            for lineno, line in enumerate(f):
                line = line.strip()
                if line.startswith(">"):
                    if header:
                        fasta_dict[header] = "".join(seq)
                        seq = []
                    header = line[1:].split(" ")[0]
                else:
                    seq.append(line)

            if header and len(seq):
                fasta_dict[header] = "".join(seq)

    except IOError as e:
        raise Exception(e)

    return fasta_dict


def write_fasta_dict(fasta_dict, filename):
    """
    Writes dictionary with fasta to file
    """
    with open(filename, "w") as f:
        for header in sorted(fasta_dict):
            f.write(">{0}\n".format(header))

            for i in range(0, len(fasta_dict[header]), 60):
                f.write(fasta_dict[header][i:i + 60] + "\n")


def parse_links_file(filename):
    scaffolds = []

    def add_contig(string):
        name = string.replace("=", "_")      #fix for nucmer
        without_sign = name[1:].strip()
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


def read_paf(filename):
    alignment = []
    for line in open(filename, "r"):
        line = line.strip()
        if not len(line):
            continue

        vals = line.split("\t")
        s_ref, e_ref = int(vals[7]), int(vals[8])
        s_qry, e_qry = int(vals[2]), int(vals[3])
        len_ref, len_qry = int(vals[6]), int(vals[1])
        ref_id, qry_id = vals[5], vals[0]

        ref_strand = 1
        qry_strand = 1
        if vals[4] == "-":
            qry_strand = -1

        ref_row = AlignmentRow(s_ref, e_ref, ref_strand, len_ref, ref_id)
        qry_row = AlignmentRow(s_qry, e_qry, qry_strand, len_qry, qry_id)
        alignment.append(AlignmentColumn(ref_row, qry_row))

    return alignment


def get_alignment(scaffolds, contigs_file, reference_file):
    query_file = "_minimap-query.fasta"
    query_fasta = {}
    contigs_fasta = read_fasta_dict(contigs_file)

    for scf in scaffolds:
        for ctg in scf.contigs:
            if "[" not in ctg.name:
                query_fasta[ctg.name] = contigs_fasta[ctg.name]
            else:
                bracket = ctg.name.index("[")
                start, end = list(map(int, ctg.name[bracket + 1: -1].split(":")))
                seq_name = ctg.name[:bracket]
                query_fasta[ctg.name] = contigs_fasta[seq_name][start:end]

    write_fasta_dict(query_fasta, query_file)

    #run minimap here
    minimap_file = "_minimap.paf"
    subprocess.check_call([MINIMAP_BIN, reference_file, query_file,  "-x", "asm5",
                           "--secondary=no", "-t", "8"], stdout=open(minimap_file, "w"))
    aln = read_paf(minimap_file)

    os.remove(query_file)
    os.remove(minimap_file)
    return aln



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


def do_job(links_file, contigs_file, reference_file):
    scaffolds = parse_links_file(links_file)
    alignment = get_alignment(scaffolds, contigs_file, reference_file)
    alignment = filter_by_coverage(alignment, 0.45)
    alignment = join_collinear(alignment)
    entry_ord, chr_len, contig_len = get_order(alignment)

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

            #flipping alignments
            for hit in entry_ord[contig.name]:
                if contig.sign < 0:
                    hit.sign = -hit.sign

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
            cur_strand = [h.sign for h in entry_ord[contig.name]]
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
            print("{0}{1}\t\t{2}\t{3}".format(sign, contig.name,
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
    if len(sys.argv) != 4:
        print("Usage: verify-minimap links contigs reference")
        return 1
    do_job(sys.argv[1], sys.argv[2], sys.argv[3])


if __name__ == "__main__":
    main()
