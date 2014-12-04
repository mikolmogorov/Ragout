
#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
Some helping function for alignment manipulations
"""

from collections import namedtuple, defaultdict
from itertools import combinations


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


##Don't use it, it's bugged!
"""
def join_collinear(alignment):
#TODO: check for strand consistency
    new_entries = []
    by_chr = group_by_chr(alignment)

    def append_entry(start_entry, last_entry):
        new_entries.append(AlignmentInfo(start_entry.ref_start,
                           last_entry.ref_end, start_entry.qry_start,
                           last_entry.qry_end,
                           abs(last_entry.ref_end - start_entry.ref_start),
                           abs(last_entry.qry_end - start_entry.qry_start),
                           last_entry.ref_id, last_entry.qry_id))

    for chr_id in by_chr:
        by_chr[chr_id].sort(key=lambda e: e.ref.start)
        start_entry = None
        last_entry = None
        for entry in by_chr[chr_id]:

            if not start_entry:
                start_entry = entry
            elif start_entry.qry_id != entry.qry_id:
                append_entry(start_entry, last_entry)
                start_entry = entry

            last_entry = entry

        if start_entry:
            append_entry(start_entry, last_entry)

    return new_entries
"""


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
        return str(self.index) + " : " + str(self.chr)


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
