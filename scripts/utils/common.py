
#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
Some helping function for alignment manipulations
"""

from __future__ import absolute_import
from collections import namedtuple, defaultdict
from six.moves import filter


AlignmentColumn = namedtuple("AlignmentColumn", ["ref", "qry"])

AlignmentRow = namedtuple("AlignmentRow", ["start", "end", "strand",
                                           "seq_len", "seq_id"])


def group_by_chr(alignment):
    by_chr = defaultdict(list)
    for entry in alignment:
        by_chr[entry.ref.seq_id].append(entry)
    for chr_id in by_chr:
        by_chr[chr_id].sort(key=lambda e: e.ref.start)
    return by_chr


def join_collinear(alignment):
    new_entries = []

    def append_entry(start_entry, last_entry):
        ref_row = AlignmentRow(start_entry.ref.start, last_entry.ref.end,
                               start_entry.ref.strand, start_entry.ref.seq_len,
                               start_entry.ref.seq_id)
        qry_left, qry_right = sorted([start_entry.qry, last_entry.qry],
                                     key=lambda r: r.start)
        qry_row = AlignmentRow(qry_left.start, qry_right.end,
                               qry_left.strand, qry_left.seq_len,
                               qry_left.seq_id)
        new_entries.append(AlignmentColumn(ref_row, qry_row))

    by_chr = group_by_chr(alignment)
    for chr_id in by_chr:
        by_chr[chr_id].sort(key=lambda e: e.ref.start)
        start_entry = None
        last_entry = None
        for entry in by_chr[chr_id]:
            if not start_entry:
                start_entry = entry
                last_entry = entry
                continue

            #checking for connection consistency
            ref_type = last_entry.ref.strand * entry.ref.strand

            qry_left, qry_right = sorted([last_entry.qry, entry.qry],
                                         key=lambda r: r.start)
            qry_type = qry_left.strand * qry_right.strand
            ###

            if (entry.qry.seq_id != last_entry.qry.seq_id or
                ref_type != qry_type):
                append_entry(start_entry, last_entry)
                start_entry = entry

            last_entry = entry

        if start_entry:
            append_entry(start_entry, last_entry)

    return new_entries


def aln_len(row):
    assert row.end >= row.start
    return row.end - row.start


def filter_by_coverage(alignment, threshold):
    by_name = defaultdict(list)
    for entry in alignment:
        by_name[entry.qry.seq_id].append(entry)

    for name in by_name:
        by_name[name].sort(key=lambda e: aln_len(e.qry), reverse=True)
        len_filter = (lambda e: aln_len(e.qry) > threshold *
                                aln_len(by_name[name][0].qry))
        by_name[name] = list(filter(len_filter, by_name[name]))

    filtered_alignment = []
    for ent_lst in by_name.values():
        filtered_alignment.extend(ent_lst)

    return filtered_alignment


class Hit:
    def __init__(self, index, chr, coord, sign):
        self.index = index
        self.chr = chr
        self.coord = coord
        self.sign = sign

    def __str__(self):
        return ("+" if self.sign > 0 else "-") + str(self.index) + " : " + str(self.chr)


def get_order(alignment):
    chr_len = defaultdict(int)
    contig_len = defaultdict(int)

    for entry in alignment:
        contig_len[entry.qry.seq_id] = max(entry.qry.end,
                                           contig_len[entry.qry.seq_id])
        chr_len[entry.ref.seq_id] = max(entry.ref.end,
                                        chr_len[entry.ref.seq_id])

    by_chr = group_by_chr(alignment)
    entry_ord = defaultdict(list)

    contig_pos = 1
    prev_start = None
    for chr_id, alignment in by_chr.items():
        for e in alignment:
            if prev_start is not None and e.ref.start > prev_start:
                contig_pos += 1
            prev_start = e.ref.start
            sign = e.ref.strand * e.qry.strand
            entry_ord[e.qry.seq_id].append(Hit(contig_pos, chr_id,
                                               e.ref.start, sign))

    return entry_ord, chr_len, contig_len
