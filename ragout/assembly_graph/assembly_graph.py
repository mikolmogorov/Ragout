#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This modules contains some functions to draw fancy graphs
(for debugging purposes mostly)
"""

import networkx as nx
import re
import logging
from collections import namedtuple

from ragout.shared import config
import assembly_refine as ar

Edge = namedtuple("Edge", ["start", "end"])
logger = logging.getLogger()

def save_colored_overlap_graph(graph_file, scaffolds, out_file):
    graph = ar._load_dot(graph_file)
    fout = open(out_file, "w")

    main_strand = set()
    all_contigs = set()
    for scf in scaffolds:
      for cont in scf.contigs:
        main_strand.add(str(cont))
        all_contigs.add(cont.name)

    fout.write("digraph {\n")
    for node in graph.nodes_iter():
        if str(node) in main_strand:
            fout.write("\"{0}\" [style=filled, fillcolor=red];\n".format(node))
        elif str(node)[1:] in all_contigs:
            fout.write("\"{0}\" [style=filled, fillcolor=blue];\n".format(node))

    for u, v in graph.edges_iter():
        fout.write("\"{0}\" -> \"{1}\";\n".format(u, v))

    fout.write("}")
    fout.close()

def save_colored_insert_overlap_graph(graph_file, scaffolds,
                                      scaffolds_refine, out_file):
    graph = ar._load_dot(graph_file)
    fout = open(out_file, "w")

    main_strand = set()
    all_contigs = set()
    for scf in scaffolds:
        main_strand |= set(map(lambda s: str(s), scf.contigs))
        all_contigs |= set(map(lambda s: s.name, scf.contigs))

    refine_contigs = set()
    for scf in scaffolds_refine:
      for cont in scf.contigs:
        if str(cont) not in main_strand:
            refine_contigs.add(str(cont))

    fout.write("digraph {\n")
    for node in graph.nodes_iter():
        if str(node) in main_strand:
            fout.write("\"{0}\" [style=filled, fillcolor=red];\n"
                       .format(node))
        elif str(node)[1:] in all_contigs:
            fout.write("\"{0}\" [style=filled, fillcolor=blue];\n"
                       .format(node))
        elif str(node) in refine_contigs:
            fout.write("\"{0}\" [style=filled, fillcolor=yellow];\n"
                       .format(node))

    for u, v in graph.edges_iter():
        fout.write("\"{0}\" -> \"{1}\";\n".format(u, v))

    fout.write("}")
    fout.close()

def save_distance_overlap_graph(graph_file, scaffolds_in, output_file):
    max_path_len = 2 * config.ASSEMBLY_MAX_PATH_LEN
    graph, edges = ar._load_dot(graph_file)

    ordered_contigs = set()
    fout = open(output_file, "w")
    fout.write("digraph {\n")
    for scf in scaffolds_in:
        for cont in scf.contigs:
            fout.write("{0} -> {1} [label=\"{2}\", color=\"red\"];\n"
                       .format(edges[str(cont)].start,
                               edges[str(cont)].end, str(cont)))
            ordered_contigs.add(str(cont))

    mark = set()
    for scf in scaffolds_in:
        for prev_cont, new_cont in zip(scf.contigs[:-1], scf.contigs[1:]):
            try:
                src = edges[str(prev_cont)].end
                dst = edges[str(new_cont)].start
            except KeyError:
                logger.debug("contigs are not in the graph")
                return None

            if src != dst and nx.has_path(graph, src, dst):
                paths = [p for p in nx.all_shortest_paths(graph, src, dst)]

                path = paths[0]
                is_good = True
                for p_start, p_end in zip(path[:-1], path[1:]):
                    found_edge = None
                    for edge_id, edge in edges.items():
                        if edge == Edge(p_start, p_end):
                            found_edge = edge_id
                            break
                    assert found_edge
                    if (found_edge in ordered_contigs and
                        found_edge != str(prev_cont) and
                        str(new_cont) != found_edge):
                        is_good = False
                        break

                if (is_good and len(path) <= 2 * max_path_len and
                    len(path) > 1 and (src, dst) not in mark):
                    mark.add((src, dst))
                    fout.write("{0} -> {1} [label=\"{2}\"];\n"
                               .format(src, dst, len(path) - 1))
    fout.write("}")
    fout.close()

def save_compress_overlap_graph(graph_file, scaffolds_in, output_file):
    graph, edges = ar._load_dot(graph_file)
    graph = nx.DiGraph(graph)

    ordered_contigs = set()
    all_contigs = set()
    for scf in scaffolds_in:
        for cont in scf.contigs:
            ordered_contigs.add(str(cont))
            all_contigs.add(cont.name)

    fout = open(output_file, "w")
    fout.write("digraph {\n")
    is_change = True
    while is_change:
        is_change = False
        for v1, v2, labels in graph.edges_iter(data=True):
            if (labels["label"] not in ordered_contigs and
                labels["label"][1:] not in all_contigs):
                is_good = True

                for y in graph.neighbors(v2):
                    if (graph[v2][y]["label"] in ordered_contigs or
                        graph[v2][y]["label"][1:] in all_contigs or y == v1):
                        is_good = False

                if is_good:
                    is_change = True
                    graph.remove_edge(v1, v2)
                    for y in graph.neighbors(v2):
                        graph.add_edge(v1, y, label=graph[v2][y]["label"])
                        graph.remove_edge(v2, y)

                if is_change:
                    break

    for v1, v2, labels in graph.edges_iter(data=True):
        if labels["label"] in ordered_contigs:
            fout.write("{0} -> {1} [label=\"{2}\", color=\"red\"];\n"
                       .format(v1, v2, labels["label"]))
        elif labels["label"][1:] in all_contigs:
            fout.write("{0} -> {1} [label=\"{2}\", color=\"blue\"];\n"
                       .format(v1, v2, labels["label"]))
        else:
            fout.write("{0} -> {1};\n".format(v1, v2))
    fout.write("}")
    fout.close()
