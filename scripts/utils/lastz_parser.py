
#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
Functions for lastz input handling
"""

from __future__ import absolute_import
from __future__ import print_function
import subprocess
import os
from itertools import combinations

from .common import AlignmentColumn, AlignmentRow
from six.moves import filter

def parse_lastz_maf(filename):
    alignments = []
    with open(filename, "r") as f:
        while True:
            line = f.readline()
            if not line:
                break

            if not line.startswith("a"):
                state = 1
                continue

            #read two next lines
            ref_line = f.readline().strip()
            qry_line = f.readline().strip()
            assert ref_line.startswith("s") and qry_line.startswith("s")

            def parse_column(string):
                seq_id, start, aln_len, strand, seq_len = string.split()[1:6]
                aln_len, seq_len = int(aln_len), int(seq_len)
                strand = 1 if strand == "+" else -1
                if strand > 0:
                    start, end = int(start), int(start) + aln_len
                else:
                    end = seq_len - 1 - int(start)
                    start = end - aln_len

                assert end >= start
                return AlignmentRow(start, end, strand, seq_len, seq_id)

            alignments.append(AlignmentColumn(parse_column(ref_line),
                                              parse_column(qry_line)))

    return alignments


def filter_intersecting(alignments):
    to_filter = set()

    def filter_by_rows(rows):
        for row_1, row_2 in combinations(rows, 2):
            if row_1.seq_id != row_2.seq_id:
                continue

            if row_1.start <= row_2.start <= row_1.end:
                to_filter.add(row_2)
                if not (row_1.start <= row_2.end <= row_1.end):
                    to_filter.add(row_1)
            continue

            #if row_2.start <= row_1.start <= row_2.end:
            #    to_filter.add(row_1)
            #    if not (row_2.start <= row_1.end <= row_2.end):
            #        to_filter.add(row_2)

    filter_by_rows([ap.ref for ap in alignments])
    filter_by_rows([ap.qry for ap in alignments])

    return [a for a in alignments if a.ref not in to_filter and
                                     a.qry not in to_filter]


def filter_by_length(alignments, min_len):
    func = (lambda a: abs(a.ref.start - a.ref.end) > min_len and
                      abs(a.qry.start - a.qry.end) > min_len)
    return list(filter(func, alignments))


LASTZ_BIN = "lastz"
def run_lastz(reference, target, out_file):
    print("Running lastz")
    reference = os.path.abspath(reference)
    target = os.path.abspath(target)
    out_file = os.path.abspath(out_file)
    cmdline = [LASTZ_BIN, reference + "[multiple,nameparse=darkspace]",
               target + "[nameparse=darkspace]","--notransition",
               "--step=20", "--chain", "--gapped", "--gfextend",
               "--ambiguous=iupac", "--format=maf", "--output=" + out_file,
               "--rdotplot=dotplot.txt"]
    devnull = os.devnull
    subprocess.check_call(cmdline, stderr=open(devnull, "w"))
