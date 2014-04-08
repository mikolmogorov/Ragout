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

    scf = Seq("")
    for line in open(ord_file, "r").read().splitlines():
        if line.startswith(">"):
            if len(scf):
                SeqIO.write(SeqRecord(scf, id=scf_name, description=""),
                            sys.stdout, "fasta")
            scf_name = line[1:]
            scf = Seq("")
        else:
            name = line[1:]
            seq = fasta_dict[name]
            if line.startswith("+"):
                scf += seq
                scf += Seq("N"  * 11)
            else:
                scf += seq.reverse_complement()
                scf += Seq("N"  * 11)
    SeqIO.write(SeqRecord(scf, id=scf_name, description=""),
                sys.stdout, "fasta")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("USAGE: ord2fasta.py ord_file contigs_file")
    else:
        main(sys.argv[1], sys.argv[2])
