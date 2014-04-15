#This module recovers an assembly graph
#by overlapping contigs
################################################

from collections import defaultdict, namedtuple
from itertools import product
import sys
import random
import logging

from Bio import SeqIO

from ragout.shared import config

logger = logging.getLogger()
Edge = namedtuple("Edge", ["begin", "end", "label"])

#PUBLIC:
#################################################

#builds assembly graph and outputs it in "dot" format
def make_overlap_graph(targets, dot_file):
    logger.info("Building overlap graph...")
    contigs_file = list(targets.values())[0]
    res = _build_overlap_graph(contigs_file, dot_file,
                               config.ASSEMBLY_MIN_OVERLAP,
                               config.ASSEMBLY_MAX_OVERLAP)
    return res


#PRIVATE:
#################################################

#reads contigs from given files
def _get_contigs(files):
    contigs = {}
    for file in files:
        for seq in SeqIO.parse(file, "fasta"):
            contigs["+" + seq.id] = str(seq.seq)
            contigs["-" + seq.id] = str(seq.seq.reverse_complement())
    return contigs

#finds overlap between two strings
def _find_overlap(str1, str2, min_ovlp, max_ovlp):
    max_k = min(max_ovlp, len(str1), len(str2))
    max_ovlp = 0
    for i in range(min_ovlp, max_k + 1):
        if str(str1)[-i:] == str(str2)[:i]:
            max_ovlp = i
    return max_ovlp


#helper function to track next unused id
def _new_node_id():
    tmp = _new_node_id.node_id
    _new_node_id.node_id += 1
    return tmp
_new_node_id.node_id = 0


def _build_overlap_graph(contigs_in, dot_out, min_ovlp, max_ovlp):
    logger.warning("C++ library for building assembly graph not found. "
                   "Using slow python version.")
    contigs = _get_contigs([contigs_in])
    edges = []
    heads = {}
    visited = set()

    dfsStack = []
    for work_ctg in contigs:
        if work_ctg in visited:
            continue

        dfsStack.append((work_ctg, _new_node_id()))
        visited.add(work_ctg)
        while len(dfsStack):
            ctg, node_id = dfsStack.pop()

            heads[ctg] = node_id
            overlaps = []
            for other_ctg in contigs:
                if ctg == other_ctg:
                    continue

                ovlp = _find_overlap(contigs[ctg], contigs[other_ctg],
                                    min_ovlp, max_ovlp)
                if ovlp >= min_ovlp:
                    overlaps.append(other_ctg)

            if overlaps:
                #ids = [heads.get(c, None) for c in overlaps]
                #if len(set(ids)) > 1:
                #    print ids

                sample_ctg = overlaps[0]
                if sample_ctg in heads:
                    cur_node = heads[sample_ctg]
                else:
                    cur_node = _new_node_id()
                    for ovlp in overlaps:
                        heads[ovlp] = cur_node
                    assert sample_ctg not in visited

                for ovlp in reversed(overlaps):
                    if not ovlp in visited:
                        dfsStack.append((ovlp, cur_node))
                        visited.add(ovlp)
            else:
                cur_node = _new_node_id()

            edges.append(Edge(node_id, cur_node, ctg))

    _output_edges(edges, dot_out)
    return True

#outputs edges to file
def _output_edges(edges, dot_file):
    fout = open(dot_file, "w")
    fout.write("digraph {\n")
    for edge in edges:
        fout.write("{0} -> {1} [label=\"{2}\"];\n".format(*edge))
    fout.write("}")

#Try load fast c++ library
try:
    from ragout.coverlap import _build_overlap_graph
except ImportError:
    pass
