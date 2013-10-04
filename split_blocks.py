#!/usr/bin/env python

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from collections import namedtuple
import sys
SyntenyBlock = namedtuple("SyntenyBlock", ["seq", "chr_id", "strand", "id", "start", "end", "chr_num"])
SeqInfo = namedtuple("SeqInfo", ["id", "length"])


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


def main(blocks_file, genomes):
    blocks, seq_info = parse_coords_file(blocks_file)
    seqs = {}
    out_files = {}
    for genome in genomes:
        seq = SeqIO.parse(genome, "fasta").next()
        seqs[seq.id] = seq.seq
        filename = genome[0:genome.rfind(".")] + "_blocks.fasta"
        out_files[seq.id] = open(filename, "w")

    for block_id, block_list in blocks.iteritems():
        if len(block_list) != len(genomes):
            continue

        for block in block_list:
            name = ""
            if block.strand == "+":
                block_seq = seqs[block.chr_id][block.start : block.end + 1]
            else:
                block_seq = seqs[block.chr_id][block.end : block.start + 1]

            #print block_id, (block.end - block.start), len(block_seq)
            block_name = "{0}_block_{1}".format(block.chr_id, block_id)
            SeqIO.write(SeqRecord(seq=block_seq, id=block_name, description=""), out_files[block.chr_id], "fasta")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
