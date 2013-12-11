from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict, namedtuple
from itertools import product
import sys
import logging

logger = logging.getLogger()

Edge = namedtuple("Edge", ["begin", "end", "label"])

def get_contigs(files):
    contigs = {}
    for file in files:
        for seq in SeqIO.parse(file, "fasta"):
            #if len(seq.seq) < 50:
            #    continue
            contigs["+" + seq.id] = seq.seq
            contigs["-" + seq.id] = seq.seq.reverse_complement()
    return contigs


#TODO: moar efficency!
def find_overlap(str1, str2, min_k):
    MAX_OVLP = 100
    max_k = min(MAX_OVLP, len(str1), len(str2))
    max_ovlp = 0
    for i in xrange(min_k, max_k + 1):
        if str(str1)[-i:] == str(str2)[:i]:
            max_ovlp = i
    return max_ovlp


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
                #for ovlp in overlaps:
                #    assert heads[ovlp] == cur_node
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


def make_overlap_graph(targets, dot_file, min_overlap):
    logger.info("Building overlap graph...")
    edges = build_graph(targets.values(), min_overlap)
    out_edges(edges, dot_file)


def out_edges(edges, dot_file):
    fout = open(dot_file, "w")
    fout.write("digraph {\n")
    for edge in edges:
        fout.write("{0} -> {1} [label=\"{2}\"];\n".format(*edge))
    fout.write("}")


if __name__ == "__main__":
    make_overlap_graph({"" : sys.argv[1]}, "overlap.dot", 21)
