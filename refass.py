#!/usr/bin/env python

import sys
import math
import graph_tools
from collections import namedtuple, defaultdict
from itertools import combinations
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


Edge = namedtuple("Edge", ["vertex", "color"])

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
    contigs = []
    permutations = []
    fin = open(filename, "r")
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
        permutations.append(map(int, blocks))

    return permutations, contigs


def build_graph(permutations):
    #find duplications
    duplications = set()
    for perm in permutations:
        current = set()
        for block in perm:
            if block in current:
                duplications.add(block)
            current.add(block)
    print duplications

    graph = defaultdict(Node)
    color = 0
    for perm in permutations:
        for i in xrange(len(perm) - 1):
            l = int(perm[i])
            r = int(perm[i + 1])
            graph[-l].edges.append(Edge(r, color))
            graph[r].edges.append(Edge(-l, color))
            color += 1
    return graph


def build_contig_index(contigs):
    index = defaultdict(list)
    for i, c in enumerate(contigs):
        for block in c.blocks:
            index[abs(block)].append(i)
    return index


def sign(val):
    return math.copysign(1, val)


def extend_scaffolds(contigs, connections):
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


def get_component_of(connected_comps, vertex):
    for con in connected_comps:
        if vertex in con:
            return con
    return None


def case_on_vs_one(graph, component, connected_comps, contig_index, num_ref):
    """
    a -- a
    """
    if len(component) != 2:
        return None

    num_edges = len(graph[component[0]].edges)
    if num_edges not in range(2, num_ref + 1):
        return None

    #print num_edges
    for fst, snd in [(0, 1), (1, 0)]:
        if abs(component[fst]) in contig_index and abs(component[snd]) not in contig_index:
            pair_comp = get_component_of(connected_comps, -component[snd])
            pair_id = pair_comp.index(-component[snd])
            other_id = abs(1 - pair_id)
            if pair_comp[other_id] in contigs:
                print "indel found!"
                return (component[fst], pair_comp[other_id])

    if abs(component[0]) in contig_index and abs(component[1]) in contig_index:
        return (component[0], component[1])

    return None


def case_indel(graph, component, connected_comps, contig_index, num_ref):
    """
    a    -b
    |  \  |
    b     c
    """
    if len(component) != 4:
        return None

    found = False
    for v1, v2 in combinations(component, 2):
        if v1 == -v2:
            found = True
            similar = [v1, v2]
            different = filter(lambda v: v != v1 and v != v2, component)

    if not found:
        return None
    #TODO: check graph structure

    #print similar, different
    if abs(similar[0]) in contig_index:
        print "deletion in some references"
        return [(s, graph[s].edges[0].vertex) for s in similar]
    else:
        print "deletion in assembly and (possibly) references"
        return [(different[0], different[1])]


def simple_connections(graph, connected_comps, contigs, num_ref):
    connections = {}
    contig_index = build_contig_index(contigs)

    for component in connected_comps:
        conn = case_on_vs_one(graph, component, conected_comps, contig_index, num_ref)
        if conn is not None:
            connections[-conn[0]] = conn[1]
            connections[-conn[1]] = conn[0]

        conn = case_indel(graph, component, conected_comps, contig_index, num_ref)
        if conn is not None:
            #print conn
            for c in conn:
                connections[-c[0]] = c[1]
                connections[-c[1]] = c[0]

    print "connections infered:", len(connections)
    return connections


def get_scaffolds(contigs, connections):
    scaffolds = extend_scaffolds(contigs, connections)
    scaffolds = filter(lambda s: len(s.contigs) > 1, scaffolds)
    for scf in scaffolds:
        for contig in scf.contigs:
            if contig.sign > 0:
                print contig.blocks,
            else:
                print map(lambda b: -b, contig.blocks)[::-1],
        print ""

    return scaffolds


def output_scaffolds(input_contigs, scaffolds, out_file):
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


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "bg.py permutations contigs"
        sys.exit(1)

    permutations, contigs = parse_permutations_file(sys.argv[1])
    graph = build_graph(permutations)
    conected_comps = graph_tools.get_connected_components(graph)
    #for c in conected_comp:
    #    print c
    connections = simple_connections(graph, conected_comps, contigs, 4)
    scaffolds = get_scaffolds(contigs, connections)

    output_scaffolds(sys.argv[2], scaffolds, "scaffolds.fasta")
    #output_graph(graph, open("bg.dot", "w"))
