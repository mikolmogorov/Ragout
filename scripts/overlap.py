#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
from collections import defaultdict
from itertools import combinations


def haming(seq1, seq2):
    count = 0
    for i in xrange(len(seq1)):
        if seq1[i] != seq2[i]:
            count += 1
    return count


def main(filename, k_size):
    kmers = defaultdict(int)
    ids = defaultdict(list)
    contigs = SeqIO.parse(filename, "fasta")
    for seq in contigs:
        km = [str(seq.seq[0:k_size]), str(seq.seq[-k_size:]),
              str(seq.seq.reverse_complement()[0:k_size]),
              str(seq.seq.reverse_complement()[-k_size:])]
        for kmer in km:
            if seq.id == "7":
                print "7", kmer
            if seq.id == "79":
                print "79", kmer
            kmers[kmer] += 1
            ids[kmer].append(seq.id)

    """
    for kmer1, kmer2 in combinations(kmers, 2):
        if haming(kmer1, kmer2) == 1:
            print "a"
    """

    count = 0
    for km in sorted(kmers):
        #print km, kmers[km]
        if kmers[km] == 1:
            assert kmers[str(Seq(km).reverse_complement())] == 1
            #print km
            print km, ids[km]
            count += 1
    print count, len(kmers)


if __name__ == "__main__":
    #main(sys.argv[1], 55)
    build_graph("final_paths.dat", "graph.dot")
