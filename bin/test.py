#!/usr/bin/env python

import sys
from collections import namedtuple, defaultdict
from itertools import product


Entry = namedtuple("Entry", ["s_ref", "e_ref", "s_qry", "e_qry",
                            "len_ref", "len_qry", "ref_id", "contig_id"])

Scaffold = namedtuple("Scaffold", ["name", "contigs"])
Contig = namedtuple("Contig", ["name", "sign", "gap"])
class Hit:
    def __init__(self, index, chr, coord):
        self.index = index
        self.chr = chr
        self.coord = coord

    def __str__(self):
        return str(self.index) + " : " + str(self.chr)


def parse_nucmer_output(filename):
    chr_alias = {}
    chr_num = 1

    entries = []
    for line in open(filename, "r"):
        line = line.strip()
        if not len(line) or not line[0].isdigit():
            continue

        vals = line.split(" | ")
        s_ref, e_ref = map(int, vals[0].split())
        s_qry, e_qry = map(int, vals[1].split())
        len_ref, len_qry = map(int, vals[2].split())
        ref_id, contig_id = vals[4].split("\t")

        if ref_id not in chr_alias:
            chr_alias[ref_id] = "chr{0}".format(chr_num)
            chr_num += 1
        entries.append(Entry(s_ref, e_ref, s_qry, e_qry,
                            len_ref, len_qry, chr_alias[ref_id], contig_id))

    return entries


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


def filter_entries(entries):
    MIN_HIT = 0.45
    by_name = defaultdict(list)
    for entry in entries:
        by_name[entry.contig_id].append(entry)

    for name in by_name:
        by_name[name].sort(key=lambda e: e.len_qry, reverse=True)
        by_name[name] = filter(lambda e: e.len_qry > MIN_HIT * by_name[name][0].len_qry, by_name[name])

    filtered_entries = []
    for ent_lst in by_name.itervalues():
        filtered_entries.extend(ent_lst)

    return filtered_entries


def join_collinear_entries(entries):
    new_entries = []
    by_chr = defaultdict(list)
    for entry in entries:
        by_chr[entry.ref_id].append(entry)
    for chr_id in by_chr:
        by_chr[chr_id].sort(key=lambda e: e.s_ref)
        #prev_contig = None
        start_entry = None
        last_entry = None
        for entry in by_chr[chr_id]:
            if not start_entry:
                start_entry = entry
            elif start_entry.contig_id != entry.contig_id:
                new_entries.append(Entry(start_entry.s_ref, last_entry.e_ref,
                                    start_entry.s_qry, last_entry.e_qry,
                                    abs(last_entry.e_ref - start_entry.s_ref),
                                    abs(last_entry.e_qry - start_entry.s_qry),
                                    last_entry.ref_id, last_entry.contig_id))
                start_entry = entry
            last_entry = entry

    return new_entries


def get_order(entries):
    chr_len = defaultdict(int)
    contig_len = defaultdict(int)

    for entry in entries:
        contig_len[entry.contig_id] = max(entry.len_qry, contig_len[entry.contig_id])
        chr_len[entry.ref_id] = max(entry.s_ref, chr_len[entry.ref_id])

    by_chr = defaultdict(list)
    for entry in entries:
        by_chr[entry.ref_id].append(entry)
    for chr_id in by_chr:
        by_chr[chr_id].sort(key=lambda e: e.s_ref)

    entry_ord = defaultdict(list)
    for chr_id, entries in by_chr.iteritems():
        for i, e in enumerate(entries):
            entry_ord[e.contig_id].append(Hit(i, chr_id, e.s_ref))

    return entry_ord, chr_len, contig_len


def gap_count(lst_1, lst_2):
    if not lst_1 or not lst_2:
        return 0

    gaps = 99999
    for i, j in product(lst_1, lst_2):
        gaps = min(gaps, abs(i.index - j.index) - 1)
    return gaps


def main():
    if len(sys.argv) < 3:
        print "Usage: test.py nucmer_coords scaffolds_ord"
        return

    entries = parse_nucmer_output(sys.argv[1])
    entries = join_collinear_entries(entries)
    entries = filter_entries(entries)
    entry_ord, chr_len, contig_len = get_order(entries)
    scaffolds = parse_contigs_order(sys.argv[2])

    #TODO: check signs
    total_breaks = 0
    total_gaps = 0
    total_contigs = 0
    for s in scaffolds:
        print "\n>" + s.name

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

            print "{0}\t{1}\t{2}".format(contig.name, contig_len[contig.name],
                                         map(str, entry_ord[contig.name]))
            total_contigs += 1

            if gap_count(prev, entry_ord[contig.name]) > 0:
                total_gaps += 1
            #total_gaps += gap_count(prev, entry_ord[contig.name])

            #only if this contig has alignments
            if entry_ord[contig.name]:
                prev = entry_ord[contig.name]

        print "\tmiss-ordered: ", len(breaks)
    print "\nTotal miss-ordered:", total_breaks
    print "Total gaps:", total_gaps
    print "Total contigs:", total_contigs
    print "Total scaffolds:", len(scaffolds)


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


if __name__ == "__main__":
    main()
