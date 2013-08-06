#!/usr/bin/env python

import sys
import math
import graph_tools
from collections import namedtuple, defaultdict
from itertools import combinations
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


Edge = namedtuple("Edge", ["vertex", "color", "distance"])
SyntenyBlock = namedtuple("SyntenyBlock", ["seq", "chr_id", "strand", "id", "start", "end", "chr_num"])
Permutation = namedtuple("Permutation", ["chr_id", "chr_num", "blocks"])
Connection = namedtuple("Connection", ["start", "end", "distance"])

class Node:
    def __init__(self):
        self.edges = []

class Contig:
    def __init__(self, name):
        self.name = name
        self.sign = 1
        self.blocks = []


class DumbContig(Contig):
    def __init__(self, length):
        Contig.__init__(self, "")
        self.length = length


class Scaffold:
    def __init__(self, left, right, contigs):
        self.left = left
        self.right = right
        self.contigs = contigs


def parse_permutations_file(filename):
    fin = open(filename, "r")
    contigs = []
    permutations = []
    contig_name = None
    ref_name = None
    ref_num = 0

    for line in fin:
        if line.startswith(">"):
            if line.startswith(">contig"):
                contig_name = line.strip()[1:]
            else:
                ref_name = line.strip()[1:]
            continue

        blocks = line.strip().split(" ")[0:-1]

        #contig
        if contig_name:
            contig = Contig(contig_name)
            contig.blocks = map(int, blocks)
            contigs.append(contig)
        #reference
        else:
            permutations.append(Permutation(chr_id=ref_name, chr_num=ref_num,
                                            blocks=map(int, blocks)))

    return (permutations, contigs)


def parse_coords_file(blocks_file):
    group = [[]]
    num_seq_id = dict()
    seq_id_num = dict()
    line = [l.strip() for l in open(blocks_file) if l.strip()]
    for l in line:
        if l[0] == '-':
            group.append([])
        else:
            group[-1].append(l)
    for l in group[0][1:]:
        l = l.split()
        num_seq_id[l[0]] = l[2]
        seq_id_num[l[2]] = int(l[0])
    ret = dict()
    for g in [g for g in group[1:] if g]:
        block_id = int(g[0].split()[1][1:])
        ret[block_id] = []
        for l in g[2:]:
            l = l.split()
            chr_id = num_seq_id[l[0]]
            start = int(l[2])
            end = int(l[3])
            chr_num = int(l[0]) - 1 #!!
            ret[block_id].append(SyntenyBlock(seq='', chr_id=chr_id, strand=l[1],
                                id=block_id, start=start, end=end, chr_num=chr_num))
    return (ret, seq_id_num)


def get_blocks_distance(left_block, right_block, ref_num, blocks_coord):
    left_instances = filter(lambda b: b.chr_num == ref_num, blocks_coord[left_block])
    right_instances = filter(lambda b: b.chr_num == ref_num, blocks_coord[right_block])

    #print len(left_instances), len(right_instances), right_instances
    assert len(left_instances) == len(right_instances) == 1
    if left_instances[0].strand == "+":
        left = left_instances[0].end
    else:
        left = left_instances[0].start

    if right_instances[0].strand == "+":
        right = right_instances[0].start
    else:
        right = right_instances[0].end

    assert right >= left
    return right - left - 1

def build_graph(permutations, blocks_coords):
    #find duplications
    duplications = set()
    for perm in permutations:
        current = set()
        for block in perm.blocks:
            if abs(block) in current:
                duplications.add(abs(block))
            current.add(abs(block))
    print duplications

    graph = defaultdict(Node)
    color = 0
    for perm in permutations:
        prev = 0
        while abs(perm.blocks[prev]) in duplications:
            prev += 1
        cur = prev + 1
        while cur < len(perm.blocks):
            while abs(perm.blocks[cur]) in duplications:
                cur += 1
            left_block = perm.blocks[prev]
            right_block = perm.blocks[cur]
            dist = get_blocks_distance(abs(left_block), abs(right_block), color, blocks_coords)
            graph[-left_block].edges.append(Edge(right_block, color, dist))
            graph[right_block].edges.append(Edge(-left_block, color, dist))
            prev = cur
            cur += 1
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


def mean(vals_list):
    assert len(vals_list) > 0
    return sum(vals_list) / len(vals_list)


def extend_scaffolds(contigs, connections):
    contig_index = build_contig_index(contigs)
    scaffolds = []
    visited = set()

    def extend_scaffold(contig):
        visited.add(contig)
        scf = Scaffold(contig.blocks[0], contig.blocks[-1], [contig])
        scaffolds.append(scf)

        #go right
        while scf.right in connections:
            adjacent = connections[scf.right].end
            distance = connections[scf.right].distance
            assert len(contig_index[abs(adjacent)]) == 1
            #    print "alarm!", len(contig_index[adjacent])
            #    break

            contig = contigs[contig_index[abs(adjacent)][0]]
            if contig in visited:
                break

            if contig.blocks[0] == adjacent:
                scf.contigs.append(DumbContig(distance))
                scf.contigs.append(contig)
                scf.right = contig.blocks[-1]
                visited.add(contig)
                continue

            if -contig.blocks[-1] == adjacent:
                scf.contigs.append(DumbContig(distance))
                scf.contigs.append(contig)
                scf.contigs[-1].sign = -1
                scf.right = -contig.blocks[0]
                visited.add(contig)
                continue

            break

        #go left
        while -scf.left in connections:
            adjacent = -connections[-scf.left].end
            distance = connections[-scf.left].distance
            assert len(contig_index[abs(adjacent)]) == 1
            #    print "alarm!", len(contig_index[adjacent])
            #    break

            contig = contigs[contig_index[abs(adjacent)][0]]
            if contig in visited:
                break

            if contig.blocks[-1] == adjacent:
                scf.contigs.insert(0, DumbContig(distance))
                scf.contigs.insert(0, contig)
                scf.left = contig.blocks[0]
                visited.add(contig)
                continue

            if -contig.blocks[0] == adjacent:
                scf.contigs.insert(0, DumbContig(distance))
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
    MIN_REF_THRESHOLD = 2   #TODO: think about it
    if len(component) != 2:
        return None

    num_edges = len(graph[component[0]].edges)
    if num_edges not in range(MIN_REF_THRESHOLD, num_ref + 1):
        return None

    #print num_edges
    for fst, snd in [(0, 1), (1, 0)]:
        if abs(component[fst]) in contig_index and abs(component[snd]) not in contig_index:
            pair_comp = get_component_of(connected_comps, -component[snd])
            pair_id = pair_comp.index(-component[snd])
            other_id = abs(1 - pair_id)
            if pair_comp[other_id] in contigs:
                print "indel found!"
                return Connection(component[fst], pair_comp[other_id], None)

    if abs(component[0]) in contig_index and abs(component[1]) in contig_index:
        start = component[0]
        end = component[1]
        distance = mean(map(lambda e:e.distance, graph[start].edges))
        return Connection(start, end, distance)

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

    if abs(similar[0]) in contig_index:
        print "deletion in some references"
        connections = []
        for s in similar:
            distance = mean(map(lambda e : e.distance, graph[s].edges))
            connections.append(Connection(s, graph[s].edges[0].vertex, distance))
        return connections
    else:
        print "deletion in assembly and (possibly) references"
        edges = filter(lambda e : e.vertex == different[1], graph[different[0]].edges)
        distance = mean(map(lambda e : e.distance, edges))
        #print distance
        return [Connection(different[0], different[1], distance)]


def simple_connections(graph, connected_comps, contigs, num_ref):
    connections = {}
    contig_index = build_contig_index(contigs)

    for component in connected_comps:
        conn = case_on_vs_one(graph, component, conected_comps, contig_index, num_ref)
        if conn is not None:
            connections[-conn.start] = Connection(-conn.start, conn.end, conn.distance)
            connections[-conn.end] = Connection(-conn.end, conn.start, conn.distance)

        conn = case_indel(graph, component, conected_comps, contig_index, num_ref)
        if conn is not None:
            for c in conn:
                connections[-c.start] = Connection(-c.start, c.end, c.distance)
                connections[-c.end] = Connection(-c.end, c.start, c.distance)

    print "connections infered:", len(connections)
    return connections


def get_scaffolds(contigs, connections):
    scaffolds = extend_scaffolds(contigs, connections)
    scaffolds = filter(lambda s: len(s.contigs) > 1, scaffolds)
    for scf in scaffolds:
        contigs = filter(lambda c : not isinstance(c, DumbContig), scf.contigs)
        for contig in contigs:
            if contig.sign > 0:
                print contig.blocks,
            else:
                print map(lambda b: -b, contig.blocks)[::-1],
        print ""

    return scaffolds


def output_scaffolds(input_contigs, scaffolds, out_file, write_contigs=False):
    contigs = SeqIO.parse(input_contigs, "fasta")
    out_stream = open(out_file, "w")
    queue = {}
    for rec in contigs:
        queue[rec.id] = rec.seq

    counter = 0
    for scf in scaffolds:
        scf_seq = Seq("")
        buffer = ""

        for i, contig in enumerate(scf.contigs):
            if isinstance(contig, DumbContig):
                buffer = "N" * contig.length
                continue

            cont_seq = queue[contig.name]
            del queue[contig.name]

            if contig.sign < 0:
                cont_seq = cont_seq.reverse_complement()

            if i > 0:
                #check for overlapping
                overlap = False
                for window in xrange(5, 100):
                    if str(scf_seq)[-window:] == str(cont_seq)[0:window]:
                        assert overlap == False
                        cont_seq = cont_seq[window:]
                        overlap = True
                if not overlap:
                    scf_seq += buffer
            buffer = ""
            scf_seq += cont_seq

        name = "scaffold{0}".format(counter)
        counter += 1
        SeqIO.write(SeqRecord(scf_seq, id=name, description=""), out_stream, "fasta")

    if write_contigs:
        for h, seq in queue.iteritems():
            SeqIO.write(SeqRecord(seq, id=h, description=""), out_stream, "fasta")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "bg.py permutations contigs"
        sys.exit(1)

    blocks_coords, seqid = parse_coords_file("data/blocks_coords.txt")
    #for r, rr in seqid.iteritems():
    #    print r, rr

    permutations, contigs = parse_permutations_file(sys.argv[1])
    num_references = len(permutations)
    graph = build_graph(permutations, blocks_coords)
    conected_comps = graph_tools.get_connected_components(graph)
    connections = simple_connections(graph, conected_comps, contigs, num_references)
    scaffolds = get_scaffolds(contigs, connections)

    output_scaffolds(sys.argv[2], scaffolds, "scaffolds.fasta")
    graph_tools.output_graph(graph, open("bg.dot", "w"))
