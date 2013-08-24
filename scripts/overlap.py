#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict, namedtuple
import sys


Node = namedtuple("Node", ["v_in", "v_out"])
def check_graph(filename, contigs):
    KMER = 55
    nodes = defaultdict(lambda : Node([], []))
    for line in open(filename, "r"):
        edge_id, start, end = line.strip().split(" ")
        nodes[start].v_out.append(edge_id)
        nodes[end].v_in.append(edge_id)

    for node in nodes.itervalues():
        kmers = []
        for edge in node.v_out:
            cont_seq = contigs[edge[1:]]
            if edge[0] == "-":
                cont_seq = cont_seq.reverse_complement()
            kmers.append(str(cont_seq)[0:KMER])
        for edge in node.v_in:
            cont_seq = contigs[edge[1:]]
            if edge[0] == "-":
                cont_seq = cont_seq.reverse_complement()
            kmers.append(str(cont_seq)[-KMER:])
        if len(kmers) > 1:
            assert len(kmers) == kmers.count(kmers[0])


def build_graph(contigs, kmer_len, out_file, out_dot):
    counter = [0]
    def new_kmer():
        counter[0] += 1
        return counter[0]

    out_dot.write("digraph {\n")
    kmer_to_num = defaultdict(new_kmer)
    for hdr, seq in contigs.iteritems():
        start = str(seq)[0:kmer_len]
        end = str(seq)[-kmer_len:]
        #out_file.write("+{0} {1} {2}\n".format(hdr, kmer_to_num[start], kmer_to_num[end]))
        out_dot.write("""{0} -> {1} [label="+{2}"]\n""".format(kmer_to_num[start], kmer_to_num[end], hdr))
        start = str(seq.reverse_complement())[0:kmer_len]
        end = str(seq.reverse_complement())[-kmer_len:]
        #out_file.write("-{0} {1} {2}\n".format(hdr, kmer_to_num[start], kmer_to_num[end]))
        out_dot.write("""{0} -> {1} [label="-{2}"]\n""".format(kmer_to_num[start], kmer_to_num[end], hdr))
    out_dot.write("}")


if __name__ == "__main__":
    contigs = {s.id : s.seq for s in SeqIO.parse(sys.argv[1], "fasta")}
    build_graph(contigs, 55, sys.stdout, open("overlap.dot", "w"))
    #check_graph(sys.argv[1], contigs)
