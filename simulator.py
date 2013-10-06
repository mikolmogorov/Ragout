#!/usr/bin/env python

import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import namedtuple
import sys
SyntenyBlock = namedtuple("SyntenyBlock", ["seq", "chr_id", "strand", "id", "start", "end", "chr_num"])
SeqInfo = namedtuple("SeqInfo", ["id", "length"])
Permutation = namedtuple("Permutation", ["chr_id", "chr_num", "blocks"])


def parse_coords_file(blocks_file):
    group = [[]]
    seq_info = {}
    blocks_info = {}
    line = [l.strip() for l in open(blocks_file) if l.strip()]
    for l in line:
        if l.startswith("-"):
            group.append([])
        else:
            group[-1].append(l)
    for l in group[0][1:]:
        seq_num, seq_len, seq_id = l.split()
        seq_num = int(seq_num)
        seq_info[seq_num] = SeqInfo(seq_id, int(seq_len))
    for g in [g for g in group[1:] if g]:
        block_id = int(g[0].split()[1][1:])
        blocks_info[block_id] = []
        for l in g[2:]:
            chr_num, bl_strand, bl_start, bl_end, bl_length = l.split()
            chr_num = int(chr_num)
            chr_id = seq_info[chr_num].id
            blocks_info[block_id].append(SyntenyBlock(seq="", chr_id=chr_id,
                                        strand=bl_strand, id=block_id,
                                        start=int(bl_start), end=int(bl_end),
                                        chr_num=int(chr_num)))
    return (blocks_info, seq_info)


def parse_permutations_file(filename, seq_info):
    fin = open(filename, "r")
    contigs = []
    permutations = []
    contig_name = None
    ref_name = None
    num_by_id = {seq.id : seq_num for seq_num, seq in seq_info.iteritems()}

    for line in fin:
        if line.startswith(">"):
            name = line.strip()[1:]
            ref_name = name
            continue

        blocks = line.strip().split(" ")[0:-1]

        ref_num = num_by_id[ref_name]
        permutations.append(Permutation(chr_id=ref_name, chr_num=ref_num,
                                        blocks=map(int, blocks)))
    return permutations


def reversal(permutation, mean_length):
    length = np.random.poisson(mean_length, 1)
    position = np.random.randint(0, len(permutation) - length)
    negate = map(lambda b: -b, permutation[position : position + length])
    return permutation[0:position] + negate[::-1] + permutation[position+length:]


def translocation(permutation, mean_length):
    length = np.random.poisson(mean_length, 1)
    pos_from = np.random.randint(0, len(permutation) - length)
    pos_to = np.random.randint(0, len(permutation) - length)

    cutted = permutation[pos_from:pos_from+length]
    rest = permutation[0:pos_from] + permutation[pos_from+length:]

    return rest[0:pos_to] + cutted + rest[pos_to:]


def evolve(permutation, branch_length):
    REVERSAL_LEN = 3
    TRANSLOC_LEN = 2
    for i in xrange(branch_length / 2):
        permutation = reversal(permutation, REVERSAL_LEN)
        permutation = translocation(permutation, TRANSLOC_LEN)
    return permutation


def main():
    print "Super Rearrangement Simulator 2000!"
    if len(sys.argv) < 4:
        print "Usage: coords_file permutations_file genome"
        return

    blocks_info, seq_info = parse_coords_file(sys.argv[1])
    permutations = parse_permutations_file(sys.argv[2], seq_info)
    sequence = SeqIO.parse(sys.argv[3], "fasta").next().seq
    root = permutations[0].blocks

    BRANCH_LEN = 6
    ref1 = evolve(root, BRANCH_LEN)
    root2 = evolve(root, BRANCH_LEN)
    ref2 = evolve(root2, BRANCH_LEN)
    root3 = evolve(root2, BRANCH_LEN)
    ref3 = evolve(root3, BRANCH_LEN)
    root4 = evolve(root3, BRANCH_LEN)
    ref4 = evolve(root4, BRANCH_LEN)
    ref5 = evolve(root4, BRANCH_LEN)

    for i, result in enumerate([ref1, ref2, ref3, ref4, ref5]):
        seq = Seq("")
        for block in result:
            binfo = filter(lambda b: b.chr_num == 1, blocks_info[abs(block)])[0]
            bseq = (sequence[binfo.start : binfo.end + 1] if binfo.strand == "+"
                    else sequence[binfo.end : binfo.start + 1])
            if block < 0:
                bseq = bseq.reverse_complement()
            seq += bseq

        ref_name = "synteinc_reference_{0}.fasta".format(i + 1)
        SeqIO.write(SeqRecord(seq=seq, id=ref_name, description=""), open(ref_name, "w"), "fasta")
    #print perm
    #print evolve(perm, 10)


if __name__ == "__main__":
    main()
