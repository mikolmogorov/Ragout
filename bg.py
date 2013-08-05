#!/usr/bin/env python

import sys
import math
from collections import namedtuple, defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

Edge = namedtuple("Edge", ["vertex", "color"])
Colors = ["red", "green", "blue", "yellow", "black"]
#Link = namedtuple("Link", ["vertex", "sign"])
#Scaffold = namedtuple("Scaffold", ["left", "right", "contigs"])


class Node:
    def __init__(self):
        self.edges = []

class Contig:
    def __init__(self, name):
        self.name = name
        self.sign = 1
        self.inverse = False
        self.blocks = []

class Scaffold:
    def __init__(self, left, right, contigs):
        self.left = left
        self.right = right
        self.contigs = contigs


def parse_permutations_file(filename):
    graph = defaultdict(Node)
    contigs = []
    fin = open(filename, "r")
    color = 0
    contig_name = None
    for line in fin:
        if line.startswith(">"):
            contig_name = line.strip()[1:] if line.startswith(">contig") else None
            continue

        blocks = line.strip().split(" ")[0:-1]
        if contig_name:
            contig = Contig(contig_name)
            contig.blocks = map(int, blocks)
            contigs.append(contig)
            continue

        #not a contig
        for i in xrange(len(blocks) - 1):
            l = int(blocks[i])
            r = int(blocks[i + 1])
            graph[-l].edges.append(Edge(r, color))
            graph[r].edges.append(Edge(-l, color))
        color += 1

    return graph, contigs


def build_contig_index(contigs):
    index = defaultdict(list)
    for i, c in enumerate(contigs):
        for block in c.blocks:
            index[abs(block)].append(i)
    return index

def sign(val):
    return math.copysign(1, val)


def get_scaffolds(contigs, connections):
    contig_index = build_contig_index(contigs)
    scaffolds = []
    visited = set()
    #for c, to in connections.iteritems():
    #    print c, to

    def extend_scaffold(contig):
        visited.add(contig)
        scf = Scaffold(contig.blocks[0], contig.blocks[-1], [contig])
        scaffolds.append(scf)

        #go right
        while scf.right in connections:
            adjacent = connections[scf.right]
            if len(contig_index[abs(adjacent)]) != 1:
                print "alarm!", len(contig_index[adjacent])
                break

            contig = contigs[contig_index[abs(adjacent)][0]]
            if contig in visited:
                break

            if contig.blocks[0] == adjacent:
                scf.contigs.append(contig)
                scf.right = contig.blocks[-1]
                visited.add(contig)
                continue

            if -contig.blocks[-1] == adjacent:
                scf.contigs.append(contig)
                scf.contigs[-1].sign = -1
                scf.right = -contig.blocks[0]
                visited.add(contig)
                continue

            break

        #go left
        while -scf.left in connections:
            adjacent = -connections[-scf.left]
            if len(contig_index[abs(adjacent)]) != 1:
                print "alarm!", len(contig_index[adjacent])
                break

            contig = contigs[contig_index[abs(adjacent)][0]]
            if contig in visited:
                break

            if contig.blocks[-1] == adjacent:
                scf.contigs.insert(0, contig)
                scf.left = contig.blocks[0]
                visited.add(contig)
                continue

            if -contig.blocks[0] == adjacent:
                scf.contigs.insert(0, contig)
                scf.contigs[0].sign = -1
                scf.left = -contig.blocks[-1]
                visited.add(contig)
                continue

            break

    for contig in contigs:
        if contig not in visited:
            extend_scaffold(contig)

    return scaffolds


def simple_connections(graph, contigs):
    NUM_REF = 4 #hardcode!
    connections = {}
    for block in graph:
        #look into breakpoint graph
        edges = map(lambda e: e.vertex, graph[block].edges)

        #only cases with connected component of 2 nodes
        if edges.count(edges[0]) == NUM_REF == len(edges):
            back_edges = map(lambda e: e.vertex, graph[edges[0]].edges)
            if back_edges.count(block) == NUM_REF == len(back_edges):
                connections[-block] = edges[0]


    scaffolds = get_scaffolds(contigs, connections)
    scaffolds = filter(lambda s: len(s.contigs) > 1, scaffolds)
    for scf in scaffolds:
        for contig in scf.contigs:
            if contig.sign > 0:
                print contig.blocks,
            else:
                print map(lambda b: -b, contig.blocks)[::-1],
        print ""

    return scaffolds


def merge_scaffolds(input_contigs, scaffolds, out_file):
    contigs = SeqIO.parse(input_contigs, "fasta")
    out_stream = open(out_file, "w")
    queue = {}
    for rec in contigs:
        queue[rec.id] = rec.seq

    counter = 0
    for scf in scaffolds:
        scf_seq = Seq("")
        for i, contig in enumerate(scf.contigs):
            cont_seq = queue[contig.name]
            del queue[contig.name]

            if contig.sign < 0:
                cont_seq = cont_seq.reverse_complement()

            if i > 0:
                #test for overlapping
                overlap = False
                for window in xrange(5, 100):
                    if str(scf_seq)[-window:] == str(cont_seq)[0:window]:
                        assert overlap == False
                        cont_seq = cont_seq[window:]
                        overlap = True
                if not overlap:
                    scf_seq += Seq("NNNNNNNNNNN")
            scf_seq += cont_seq

        name = "scaffold{0}".format(counter)
        counter += 1
        SeqIO.write(SeqRecord(scf_seq, id=name, description=""), out_stream, "fasta")

    #for h, seq in queue.iteritems():
    #    SeqIO.write(SeqRecord(seq, id=h, description=""), out_stream, "fasta")


def output_graph(graph, dot_file):
    dot_file.write("graph {\n")
    used_vertexes = set()
    for node_id, node in graph.iteritems():
        for edge in node.edges:
            if edge.vertex not in used_vertexes:
                dot_file.write("""{0} -- {1} [color = "{2}"];\n"""
                                .format(node_id, edge.vertex, Colors[edge.color]))
        used_vertexes.add(node_id)

    #for i in xrange(max(graph.keys()) + 1):
    #    dot_file.write("""{0} -- {1} [color = "black"];\n""".format(i, -i))
    dot_file.write("}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "bg.py permutations contigs"
        sys.exit(1)

    graph, contigs = parse_permutations_file(sys.argv[1])
    scaffolds = simple_connections(graph, contigs)
    merge_scaffolds(sys.argv[2], scaffolds, "scaffolds.fasta")
    #output_graph(graph, open("bg.dot", "w"))
