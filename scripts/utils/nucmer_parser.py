
from .common import AlignmentRow, AlignmentColumn
from six.moves import map

#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

#Some helper functions to parse/process nucmer output

def parse_nucmer_coords(filename):
    chr_alias = {}
    chr_num = 1

    alignment = []
    for line in open(filename, "r"):
        line = line.strip()
        if not len(line) or not line[0].isdigit():
            continue

        vals = line.split(" | ")
        s_ref, e_ref = list(map(int, vals[0].split()))
        s_qry, e_qry = list(map(int, vals[1].split()))
        len_ref, len_qry = list(map(int, vals[2].split()))
        ref_id, qry_id = vals[4].split("\t")

        if e_ref > s_ref:
            ref_strand = 1
        else:
            ref_strand = -1
            s_ref, e_ref = e_ref, s_ref

        if e_qry > s_qry:
            qry_strand = 1
        else:
            qry_strand = -1
            s_qry, e_qry = e_qry, s_qry

        if ref_id not in chr_alias:
            chr_alias[ref_id] = "chr{0}".format(chr_num)
            chr_num += 1

        ref_row = AlignmentRow(s_ref, e_ref, ref_strand, None,
                                                chr_alias[ref_id])
        qry_row = AlignmentRow(s_qry, e_qry, qry_strand, None, qry_id)
        alignment.append(AlignmentColumn(ref_row, qry_row))

    return alignment
