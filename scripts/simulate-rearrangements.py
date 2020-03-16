#!/usr/bin/env python

#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
A script that simulates inversions in a given genome.
"""

from __future__ import absolute_import
from __future__ import print_function
import sys
import random
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from utils.nucmer_parser import *
from six.moves import range


def get_unique_contigs(alignment):
    first_filtration = defaultdict(int)
    for e in alignment:
        first_filtration[e.contig_id] += 1

    result = []
    for cont, counts in first_filtration.items():
        if counts == 1:
            result.append(cont)

    return result


def get_contigs_with_length(unique_names, alignment_coords, treshold):
    result = []
    for name in unique_names:
        start, end  = alignment_coords[name]
        if abs(start - end) >= treshold:
            result.append(name)
    return result


def do_job(nucmer_coords, number_of_inv, orig_reference,
           output_reference, treshold):
    alignment = parse_nucmer_coords(nucmer_coords)
    alignment = join_collinear(alignment)
    alignment = filter_by_coverage(alignment)

    unique_seq = get_unique_contigs(alignment)

    alignment_coords = {}
    for e in alignment:
        alignment_coords[e.contig_id] = (int(e.s_ref), int(e.e_ref))

    unique_seq = get_contigs_with_length(unique_seq, alignment_coords,
                                         treshold)

    seq = list(SeqIO.parse(orig_reference, "fasta"))[0]
    refer_name = seq.name
    refer = seq.seq

    for _ in range(number_of_inv):
        contig = random.choice(unique_seq)
        start, end = alignment_coords[contig]
        compl = refer[start : end].reverse_complement()
        refer = refer[:start] + compl + refer[end:]

        print((contig + " (" + str(start) + ", " + str(end) +  ")"))

        unique_seq.remove(contig)
        if not unique_seq:
            break

    SeqIO.write(SeqRecord(id=refer_name, description="", seq=refer),
                output_reference, "fasta")


def main():
    if len(sys.argv) != 6:
        print("Usage: simulate-rearrangements.py nucmer_coords "
              "orig_ref inv_num min_length out_ref\n\n"
              "A script which creates a new reference with simulated "
              "inversions from a given one. Made for testing purposes. "
              "This script requires an original reference and a nucmer alignment "
              "of contigs (in coords format) on that reference. One should use "
              "contigs whcich will be used to run Ragout in future. "
              "Reference is expected to have "
              "only one fasta sequence (one chromosome).\n\n"
              "Positional arguments:\n"
              "nucmer_coords\tcontigs alignment on original reference\n"
              "orig_ref\tpath to original reference (will be transformed)\n"
              "inv_num\t\tnumber of inversions to simulate\n"
              "min_length\tminimum length of inversion\n"
              "out_ref\t\tpath to reference that will be created")
        return

    nucmer_coords = sys.argv[1]
    orig_reference = sys.argv[2]
    number_of_inv = int(sys.argv[3])
    treshold = int(sys.argv[4])
    output_reference = sys.argv[5]

    do_job(nucmer_coords, number_of_inv, orig_reference,
           output_reference, treshold)

if __name__ == "__main__":
    main()
