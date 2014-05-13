#!/usr/bin/env python

from __future__ import print_function
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys

def main(ord_file, contigs_file):
    fasta_dict = {}
    for seq in SeqIO.parse(contigs_file, "fasta"):
        fasta_dict[seq.id] = seq.seq

    scf_seqs = []
    for line in open(ord_file, "r").read().splitlines():
        if line.startswith(">"):
            if len(scf_seqs):
                scf_str = "".join(scf_seqs)
                SeqIO.write(SeqRecord(Seq(scf_str), id=scf_name,
                            description=""), sys.stdout, "fasta")
            scf_name = line[1:]
            scf_seqs = []
        else:
            name = line[1:]
            seq = fasta_dict[name]
            if line.startswith("+"):
                scf_seqs.append(str(seq))
                scf_seqs.append("N" * 11)
            else:
                scf_seqs.append(str(seq.reverse_complement()))
                scf_seqs.append("N" * 11)
    scf_str = "".join(scf_seqs)
    SeqIO.write(SeqRecord(Seq(scf_str), id=scf_name, description=""),
                sys.stdout, "fasta")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("USAGE: ord2fasta.py ord_file contigs_file")
    else:
        main(sys.argv[1], sys.argv[2])
