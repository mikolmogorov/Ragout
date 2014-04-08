#!/usr/bin/env python

from __future__ import print_function
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os, sys
import subprocess
from collections import namedtuple

Block = namedtuple("Block", ["id", "start", "length", "chr_id"])
SeqInfo = namedtuple("SeqInfo", ["id", "length"])


GEPARD_DIR = "/home/volrath/Bioinf/tools/gepard-1.30/"
EXEC = "./gepardcmd.sh"
MATRIX = "matrices/edna.mat"

def main():
    if len(sys.argv) < 4:
        print("Usage: dotplot.py <coords_file> <output_dir> <fasta_references>")
        return
    coords_file = sys.argv[1]
    fasta = sys.argv[3:]
    blocks = parse_coords_file(coords_file)
    draw_dot_plot(blocks, fasta, sys.argv[2])


def draw_dot_plot(blocks, seq_files, out_dir):
    seqs = {}
    for file in seq_files:
        for seq in SeqIO.parse(file, "fasta"):
            seq_id = seq.id.split(" ")[0]
            seqs[seq_id] = seq.seq

    out_dir = os.path.abspath(out_dir)
    os.chdir(GEPARD_DIR)

    for block_id, blocklist in blocks.items():
        allBlocks = open(os.path.join(out_dir, "blocks{0}.fasta".format(block_id)), "w")
        for block in blocklist:
            seq = seqs[block.chr_id][block.start : block.start + block.length]
            if block.id < 0:
                seq = seq.reverse_complement()
            SeqIO.write(SeqRecord(seq, id=block.chr_id, description=""),
                        allBlocks, "fasta")


        block1, block2 = blocklist[0:2]

        s1 = seqs[block1.chr_id]
        s2 = seqs[block2.chr_id]

        seq1 = s1[block1.start : block1.start + block1.length]
        seq2 = s2[block2.start : block2.start + block2.length]
        if block1.id < 0:
            seq1 = seq1.reverse_complement()
        if block2.id < 0:
            seq2 = seq2.reverse_complement()


        file1 = os.path.join(out_dir, "block{0}.{1}-1.fasta".format(block_id, block1.chr_id))
        file2 = os.path.join(out_dir, "block{0}.{1}-2.fasta".format(block_id, block2.chr_id))

        SeqIO.write(SeqRecord(seq1), file1, "fasta")
        SeqIO.write(SeqRecord(seq2), file2, "fasta")

        out_dot = os.path.join(out_dir, "block{0}-{1}-{2}.png"
                                            .format(block_id, block1.chr_id, block2.chr_id))

        cmdline = [EXEC, "-seq1", file1, "-seq2", file2, "-outfile",
                    out_dot, "-matrix", MATRIX, "-word", "10"]
        subprocess.check_call(cmdline)

        os.remove(file1)
        os.remove(file2)


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

            num_id = block_id if bl_strand == "+" else -block_id
            blocks_info[block_id].append(Block(num_id, int(bl_start),
                                        int(bl_length), chr_id))
    return blocks_info

if __name__ == "__main__":
    main()
