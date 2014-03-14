#This module recovers an assembly graph
#by overlapping contigs
################################################

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict, namedtuple
from itertools import product
import sys
import random
import logging

import config

logger = logging.getLogger()
Edge = namedtuple("Edge", ["begin", "end", "label"])

#PUBLIC:
#################################################

#builds assembly graph and outputs it in "dot" format
def make_overlap_graph(targets, dot_file):
    logger.info("Building overlap graph...")
    #edges = build_graph(targets.values(), config.ASSEMBLY_MIN_OVERLAP)
    contigs = get_contigs(targets.values())
    kmer = guess_kmer(contigs)
    logger.info("Kmer value seems to be equal to " + str(kmer))
    edges = build_with_hash(contigs, kmer)
    out_edges(edges, dot_file)


#PRIVATE:
#################################################

#reads contigs from given files
def get_contigs(files):
    contigs = {}
    for file in files:
        for seq in SeqIO.parse(file, "fasta"):
            contigs["+" + seq.id] = str(seq.seq)
            contigs["-" + seq.id] = str(seq.seq.reverse_complement())
    return contigs


def build_with_hash(contigs, kmer_len):
    next_kmer_id = [0]
    kmers_ids = {}
    def kmer_enum(kmer):
        if not kmer in kmers_ids:
            kmers_ids[kmer] = next_kmer_id[0]
            next_kmer_id[0] += 1
        return kmers_ids[kmer]

    edges = []
    for hdr, seq in contigs.iteritems():
        begin = kmer_enum(seq[:kmer_len])
        end = kmer_enum(seq[-kmer_len:])
        edges.append(Edge(begin, end, hdr))

    return edges


#finds overlap between two strings
def find_overlap(str1, str2, min_k):
    max_k = min(config.ASSEMBLY_MAX_OVERLAP, len(str1), len(str2))
    max_ovlp = 0
    for i in xrange(min_k, max_k + 1):
        if str(str1)[-i:] == str(str2)[:i]:
            max_ovlp = i
    return max_ovlp


def most_common(lst):
    return max(set(lst), key=lst.count)


def guess_kmer(contigs):
    to_compare = random.sample(contigs.keys(), 100)
    overlaps = []
    for hdr1 in to_compare:
        for hdr2 in contigs:
            if hdr1 == hdr2:
                continue
            ovlp = find_overlap(contigs[hdr1], contigs[hdr2],
                                10)
            if ovlp:
                overlaps.append(ovlp)
    print overlaps
    return most_common(overlaps)


"""
#helper function to track next unused id
def new_node_id():
    tmp = new_node_id.node_id
    new_node_id.node_id += 1
    return tmp
new_node_id.node_id = 0


def build_graph(files, min_ovlp):
    contigs = get_contigs(files)
    edges = []
    heads = {}
    visited = set()

    def dfs(ctg, node_id):
        visited.add(ctg)
        heads[ctg] = node_id
        overlaps = []
        for other_ctg in contigs:
            if ctg == other_ctg:
                continue

            ovlp = find_overlap(contigs[ctg], contigs[other_ctg], min_ovlp)
            if ovlp:
                overlaps.append(other_ctg)

        if overlaps:
            sample_ctg = overlaps[0]
            if sample_ctg in heads:
                cur_node = heads[sample_ctg]
            else:
                cur_node = new_node_id()
                for ovlp in overlaps:
                    heads[ovlp] = cur_node
                assert sample_ctg not in visited

            for ovlp in overlaps:
                if not ovlp in visited:
                    dfs(ovlp, cur_node)
        else:
            cur_node = new_node_id()

        edges.append(Edge(node_id, cur_node, ctg))

    for ctg in contigs:
        if not ctg in visited:
            dfs(ctg, new_node_id())

    return edges
"""

#outputs edges to file
def out_edges(edges, dot_file):
    fout = open(dot_file, "w")
    fout.write("digraph {\n")
    for edge in edges:
        fout.write("{0} -> {1} [label=\"{2}\"];\n".format(*edge))
    fout.write("}")
