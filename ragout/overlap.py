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
    edges = build_graph(targets.values(), config.ASSEMBLY_MIN_OVERLAP)
    #contigs = get_contigs(targets.values())
    #kmer = 63#guess_kmer(contigs)
    #logger.info("Kmer value seems to be equal to " + str(kmer))
    #edges = build_with_hash(contigs, kmer)
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
    n_overlap = [0]
    def kmer_enum(kmer):
        if not kmer in kmers_ids:
            kmers_ids[kmer] = next_kmer_id[0]
            next_kmer_id[0] += 1
        else:
            n_overlap[0] += 1
        return kmers_ids[kmer]

    edges = []
    for hdr, seq in contigs.iteritems():
        begin = kmer_enum(seq[:kmer_len])
        end = kmer_enum(seq[-kmer_len:])
        edges.append(Edge(begin, end, hdr))

    print n_overlap[0]
    return edges

#############
def prefix_func(string):
    func = [0] * len(string)
    for i in xrange(1, len(string)):
        k = func[i - 1]
        while True:
            if string[i] == string[k]:
                func[i] = k + 1
                break
            if k == 0:
                func[i] = 0
                break
            k = func[k - 1]
    return func


def kmp_overlap(string_head, string_tail):
    pfun = prefix_func(string_tail)
    i, j = 0, 0
    #found = []
    max_ovlp = 0
    while i + j < len(string_head):
        if string_head[i + j] == string_tail[j]:
            if i + j == len(string_head) - 1:
                #print "found overlap by", j + 1
                max_ovlp = max(j + 1, max_ovlp)

            if j == len(string_tail) - 1:
                #found.append(i)
                i = i + j - pfun[j - 1]
                j = pfun[j - 1]
            else:
                j += 1

        else:
            if j == 0:
                i += 1
            else:
                i = i + j - pfun[j - 1]
                j = pfun[j - 1]

    return max_ovlp
######################


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
            ovlp = find_overlap(contigs[hdr1], contigs[hdr2], 10)
            if ovlp:
                overlaps.append(ovlp)
    print overlaps
    return most_common(overlaps)


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
    #n_overlaps = [0]
    #counter = 0

    dfsStack = []
    for work_ctg in contigs:
        if work_ctg in visited:
            continue

        dfsStack.append((work_ctg, new_node_id()))
        visited.add(work_ctg)
        while len(dfsStack):
            ctg, node_id = dfsStack.pop()

            heads[ctg] = node_id
            overlaps = []
            #counter += 1
            for other_ctg in contigs:
                if ctg == other_ctg:
                    continue

                #ovlp = find_overlap(contigs[ctg], contigs[other_ctg], min_ovlp)
                max_k = config.ASSEMBLY_MAX_OVERLAP
                ovlp = kmp_overlap(contigs[ctg][-max_k:], contigs[other_ctg][:max_k])
                if ovlp > min_ovlp:
                    overlaps.append(other_ctg)
                    #n_overlaps[0] += 1

            if overlaps:
                sample_ctg = overlaps[0]
                if sample_ctg in heads:
                    cur_node = heads[sample_ctg]
                else:
                    cur_node = new_node_id()
                    for ovlp in overlaps:
                        heads[ovlp] = cur_node
                    assert sample_ctg not in visited

                for ovlp in reversed(overlaps):
                    if not ovlp in visited:
                        dfsStack.append((ovlp, cur_node))
                        visited.add(ovlp)
            else:
                cur_node = new_node_id()

            edges.append(Edge(node_id, cur_node, ctg))


    #print n_overlaps[0]
    #print counter
    return edges


#outputs edges to file
def out_edges(edges, dot_file):
    fout = open(dot_file, "w")
    fout.write("digraph {\n")
    for edge in edges:
        fout.write("{0} -> {1} [label=\"{2}\"];\n".format(*edge))
    fout.write("}")
