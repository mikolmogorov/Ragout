import networkx as nx
import sys
import re
from collections import namedtuple
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from datatypes import Contig, Scaffold


Edge = namedtuple("Edge", ["start", "end"])


#in order not to depend on heavy graphviz
def load_dot(filename):
    graph = nx.DiGraph()
    edges = {}
    pattern = re.compile("([0-9]+)\s*\->\s*([0-9]+)\s*\[.*=.*\"(.+)\".*\];")
    for line in open(filename, "r").read().splitlines():
        m = pattern.match(line)
        if not m:
            continue

        v1, v2 = int(m.group(1)), int(m.group(2))
        graph.add_node(v1)
        graph.add_node(v2)
        graph.add_edge(v1, v2, label=m.group(3))
        edges[m.group(3)] = Edge(v1, v2)
    return graph, edges


def insert_from_graph(graph_file, scaffolds_in):
    MAX_LEN = 6
    new_scaffolds = []
    graph, edges = load_dot(graph_file)

    ordered_contigs = set()
    for scf in scaffolds_in:
        ordered_contigs |= set(map(lambda s: s.name, scf.contigs))

    for scf in scaffolds_in:
        new_scaffolds.append(Scaffold(scf.name))

        for prev_cont, new_cont in zip(scf.contigs[:-1], scf.contigs[1:]):
            new_scaffolds[-1].contigs.append(prev_cont)

            #print str(prev_cont)
            try:
                src = edges[str(prev_cont)].end
                dst = edges[str(new_cont)].start
            except KeyError:
                #print "contig is not in graph"
                continue

            if src == dst:
                #print "adjacent contigs {0} and {1}".format(prev_cont, new_cont)
                new_scaffolds[-1].contigs[-1].gap = 0
                continue

            if not nx.has_path(graph, src, dst):
                #print "no path between {0} and {1}".format(prev_cont, new_cont)
                continue

            paths = [p for p in nx.all_shortest_paths(graph, src, dst)]
            if len(paths) != 1:
                #print "ambigious paths between {0} and {1}".format(prev_cont, new_cont)
                continue

            path = paths[0]
            #print len(path)
            if len(path) > MAX_LEN:
                #print "too long path between {0} and {1}".format(prev_cont, new_cont)
                continue

            #all is ok!
            print ("found path between {0} and {1} of length {2}"
                                    .format(prev_cont, new_cont, len(path)))
            for p_start, p_end in zip(path[:-1], path[1:]):
                #corresponging edge in graph
                found_edge = None
                for edge_id, edge in edges.iteritems():
                    if edge == Edge(p_start, p_end):
                        found_edge = edge_id
                        break
                assert found_edge

                if found_edge[1:] in ordered_contigs:
                    print "Alarm! path inconsistency:", found_edge

                new_scaffolds[-1].contigs[-1].gap = 0
                new_scaffolds[-1].contigs.append(Contig.from_sting(found_edge))
                new_scaffolds[-1].contigs[-1].gap = 0

        new_scaffolds[-1].contigs.append(new_cont)

    return new_scaffolds


def refine_contigs(graph_file, scaffolds):
    new_scaffolds = insert_from_graph(graph_file, scaffolds)
    return new_scaffolds
