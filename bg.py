#!/usr/bin/env python

import sys
import math
from collections import namedtuple, defaultdict

Edge = namedtuple("Edge", ["vertex", "color"])
Colors = ["red", "green", "blue", "yellow", "black"]
Link = namedtuple("Link", ["vertex", "sign"])
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

    def extend_scaffold(contig):
        visited.add(contig)
        scf = Scaffold(contig.blocks[0], contig.blocks[-1], [contig])
        scaffolds.append(scf)

        #go right
        while abs(scf.right) in connections:
            adjacent = connections[abs(scf.right)]
            assert len(contig_index[adjacent.vertex]) == 1
            contig = contigs[contig_index[adjacent.vertex][0]]
            if contig in visited:
                break

            if abs(contig.blocks[0]) == adjacent.vertex:
                scf.contigs.append(contig)
                scf.contigs[-1].sign = sign(scf.right * contig.blocks[0] * adjacent.sign)
                scf.right = contig.blocks[-1] * scf.contigs[-1].sign
                visited.add(contig)
            elif abs(contig.blocks[-1]) == adjacent.vertex:
                scf.contigs.append(contig)
                scf.contigs[-1].sign = sign(scf.right * contig.blocks[-1] * adjacent.sign)
                scf.contigs[-1].inverse = True
                scf.right = contig.blocks[0] * scf.contigs[-1].sign
                visited.add(contig)

        #go left
        while abs(scf.left) in connections:
            adjacent = connections[abs(scf.left)]
            assert len(contig_index[adjacent.vertex]) == 1
            contig = contigs[contig_index[adjacent.vertex][0]]
            if contig in visited:
                break

            if abs(contig.blocks[-1]) == adjacent.vertex:
                scf.contigs.insert(0, contig)
                scf.contigs[0].sign = sign(scf.left * contig.blocks[-1] * adjacent.sign)
                scf.left = contig.blocks[0] * scf.contigs[0].sign
                visited.add(contig)
            elif abs(contig.blocks[0]) == adjacent.vertex:
                scf.contigs.insert(0, contig)
                scf.contigs[0].sign = sign(scf.left * contig.blocks[0] * adjacent.sign)
                scf.contigs[0].inverse = True
                scf.left = contig.blocks[-1] * scf.contigs[0].sign
                visited.add(contig)


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
                s = sign(-block * edges[0])
                connections[abs(block)] = Link(vertex=abs(edges[0]), sign=s)
                connections[abs(edges[0])] = Link(vertex=abs(block), sign=s)


    for scf in get_scaffolds(contigs, connections):
        for contig in scf.contigs:
            if not contig.inverse:
                print map(lambda b: int(b * contig.sign), contig.blocks),
            else:
                print map(lambda b: int(b * contig.sign), contig.blocks)[::-1],
        print ""


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
    graph, contigs = parse_permutations_file(sys.argv[1])
    simple_connections(graph, contigs)
    #output_graph(graph, open("bg.dot", "w"))
