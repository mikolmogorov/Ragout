
import networkx as nx
import re
import logging
from collections import namedtuple

import assembly_refine as ar

Edge = namedtuple("Edge", ["start", "end"])
logger = logging.getLogger()

def save_colored_overlap_graph(graph_file, scaffolds, out_file):
    graph, edges = ar._load_dot(graph_file)
    fout = open(out_file, "w")

    main_strand = set()
    all_contigs = set()
    for scf in scaffolds:
      for cont in scf.contigs:
        main_strand.add(str(cont))
        all_contigs.add(cont.name)

    fout.write("digraph {\n")
    for edge_id, edge in edges.iteritems():
        if edge_id in main_strand:
            fout.write("{0} -> {1} [label=\"{2}\", color=\"red\"];\n".format((edge_id, edge.start, edge.end)))
        elif edge_id[1:] in all_contigs:
            fout.write("{0} -> {1} [label=\"{2}\", color=\"blue\"];\n".format((edge_id, edge.start, edge.end)))
        else:
            fout.write("{0} -> {1} [label=\"{2}\"];\n".format((edge_id, edge.start, edge.end)))
    fout.write("}")

def save_compress_overlap_graph(graph_file, scaffolds, output_file):
    pass

'''
def compress_overlap_graph(graph_file, scaffolds, output_file):
    max_path_len = 2 * config.ASSEMBLY_MAX_PATH_LEN
    graph, edges = ar._load_dot(graph_file)
    compress_edges = []

    ordered_contigs = set()
    for scf in scaffolds_in:
        ordered_contigs |= set(map(lambda s: s.name, scf.contigs))

    for scf in scaffolds_in:
        new_scaffolds.append(Scaffold(scf.name))

        for prev_cont, new_cont in zip(scf.contigs[:-1], scf.contigs[1:]):
            try:
                src = edges[str(prev_cont)].end
                dst = edges[str(new_cont)].start
            except KeyError:
                logger.debug("contigs are not in the graph")
                return None

        if src != dst and nx.has_path(graph, src, dst):
            path_edges = [p for p in nx.all_shortest_paths(graph, src, dst)]
'''