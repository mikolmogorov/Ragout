#!/usr/bin/env python

from collections import namedtuple
import sys
import scripts.debrujin_refine as debrujin
from scripts.datatypes import *


def parse_contigs_order(filename):
    scaffolds = []
    for line in open(filename, "r"):
        if line.startswith(">"):
            scaffolds.append(Scaffold(line.strip()[1:]))
        else:
            name = line.strip("\n").replace("=", "_") #fix for quast
            without_sign = name[1:].strip()
            sign = 1 if name[0] == "+" else -1
            scaffolds[-1].contigs.append(Contig(without_sign, sign))
    return scaffolds


def merge(big_file, small_file):
    big_scaffolds = parse_contigs_order(big_file)
    small_scaffolds = parse_contigs_order(small_file)

    big_index = set()
    for scf in big_scaffolds:
        for c in scf.contigs:
            big_index.add(c.name)

    count = 0
    new_scafflods = []
    for scf in big_scaffolds:
        result = []
        for prev_cont, new_cont in zip(scf.contigs[:-1], scf.contigs[1:]):
            found_pair = False
            for small_scf in small_scaffolds:
                names = map(lambda c: c.name, small_scf.contigs)
                try:
                    begin = names.index(prev_cont.name)
                    end = names.index(new_cont.name)
                    found_pair = True
                    break
                except ValueError:
                    continue

            result.append(prev_cont)
            if not found_pair:
                continue

            assert end != begin
            same_dir = True
            if end < begin:
                same_dir = False
                end, begin = begin, end

            consistent = True
            for c in small_scf.contigs[begin + 1 : end]:
                if c.name in big_index:
                    consistent = False
                    break

            if not consistent or end - begin == 1:
                continue

            if ((prev_cont.sign == new_cont.sign) !=
                (small_scf.contigs[begin].sign == small_scf.contigs[end].sign)):
                continue

            count += end - begin - 1
            contigs = small_scf.contigs[begin + 1 : end]
            if not same_dir:
                contigs = contigs[::-1]
                contigs = map(lambda c: Contig(c.name, -c.sign, 0), contigs)
            result.extend(contigs)

        result.append(new_cont)
        s = Scaffold(scf.name)
        s.contigs = result
        new_scafflods.append(s)

    return new_scafflods

def main():
    if len(sys.argv) < 4:
        print "Usage: merge_iters.py big small overlap"

    scaffolds = merge(sys.argv[1], sys.argv[2])
    refined_scaffolds = debrujin.refine_contigs(sys.argv[3], scaffolds)

    scaf_out = open("merged.txt", "w")
    ref_out = open("merged_refined.txt", "w")

    for scf in scaffolds:
        scaf_out.write(">" + scf.name + "\n")
        for c in scf.contigs:
            sign = "+" if c.sign > 0 else "-"
            scaf_out.write(sign + c.name + "\n")

    for scf in refined_scaffolds:
        ref_out.write(">" + scf.name + "\n")
        for c in scf.contigs:
            sign = "+" if c.sign > 0 else "-"
            ref_out.write(sign + c.name + "\n")


if __name__ == "__main__":
    main()
