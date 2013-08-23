#!/usr/bin/end python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys


def do(path_file, contigs_file, out_path, out_contigs):
    contigs = SeqIO.parse(contigs_file, "fasta")

    p_out = open(out_path, "w")
    c_out = open(out_contigs, "w")

    counter = 1
    conj = False
    for line in open(path_file, "r"):
        if not line.startswith("PATH"):
            continue

        vals = line.strip().split(" ")
        edge_id, conj_id, start, end = vals[1], vals[6], int(vals[8]), int(vals[10])

        if not conj:
            p_out.write("+seq{0} {1} {2}\n".format(counter, start, end))
            record = contigs.next()
            SeqIO.write(SeqRecord(id="seq{0}".format(counter),
                        seq=record.seq, description=""), c_out, "fasta")
            conj = True
        else:
            record = contigs.next()
            p_out.write("-seq{0} {1} {2}\n".format(counter, start, end))
            counter += 1
            conj = False

    assert not conj


def main():
    if len(sys.argv) < 3:
        print "Usage: preprocess.py path_file contigs_file"
        return
    do(sys.argv[1], sys.argv[2], "path.txt", "path.fasta")


if __name__ == "__main__":
    main()
