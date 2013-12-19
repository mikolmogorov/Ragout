#!/usr/bin/env python

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from collections import namedtuple
import sys, os
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


def check_block(ngenomes, block_list):
    MIN = 2000
    MAX = 5000
    #ids = map(lambda x: x.chr_id, block_list)
    if check_unique(block_list) and len(block_list) == ngenomes:
        length = abs(block_list[0].start - block_list[0].end)
        if MIN <= length <= MAX:
            return True
    return False


def check_unique(block_list):
    ids = map(lambda x: x.chr_id, block_list)
    return len(ids) == len(set(ids))


def main(blocks_file, genomes):
    blocks, seq_info = parse_coords_file(blocks_file)
    seqs = {}
    out_files = {}
    id_to_name = {}
    for genome in genomes:
        seq = SeqIO.parse(genome, "fasta").next()
        seqs[seq.id] = seq.seq
        filename = genome[0:genome.rfind(".")] + "_blocks.fasta"
        out_files[seq.id] = open(filename, "w")
        id_to_name[seq.id] = os.path.basename(genome)

    for block_id, block_list in blocks.iteritems():
        #if len(block_list) != len(genomes):
        #    continue
        if not check_unique(block_list):
            continue

        """
        if check_block(len(genomes), block_list):
            out_block = open("block_" + str(block_id) + ".fasta", "w")
            for block in block_list:
                if block.strand == "+":
                    block_seq = seqs[block.chr_id][block.start : block.end + 1]
                else:
                    block_seq = seqs[block.chr_id][block.end : block.start + 1]

                name = id_to_name[block.chr_id]
                SeqIO.write(SeqRecord(seq=block_seq, id=name, description=""), out_block, "fasta")
        """

        for block in block_list:
            name = ""
            if block.strand == "+":
                block_seq = seqs[block.chr_id][block.start : block.end + 1]
            else:
                block_seq = seqs[block.chr_id][block.end : block.start + 1]

            block_name = "{0}_block_{1}".format(block.chr_id, block_id)
            SeqIO.write(SeqRecord(seq=block_seq, id=block_name, description=""), out_files[block.chr_id], "fasta")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
