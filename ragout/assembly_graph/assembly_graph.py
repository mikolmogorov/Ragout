
import networkx as nx
import re
import logging
from collections import namedtuple

from shared import config
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
            fout.write("{0} -> {1} [label=\"{2}\", color=\"red\"];\n".format(edge.start, edge.end, edge_id))
        elif edge_id[1:] in all_contigs:
            fout.write("{0} -> {1} [label=\"{2}\", color=\"blue\"];\n".format(edge.start, edge.end, edge_id))
        else:
            fout.write("{0} -> {1} [label=\"{2}\"];\n".format(edge.start, edge.end, edge_id))
    fout.write("}")

def save_colored_insert_overlap_graph(graph_file, scaffolds, scaffolds_refine, out_file):
    graph, edges = ar._load_dot(graph_file)
    fout = open(out_file, "w")

    main_strand = set()
    all_contigs = set()
    for scf in scaffolds:
      for cont in scf.contigs:
        main_strand.add(str(cont))
        all_contigs.add(cont.name)

    refine_contigs = set()
    for scf in scaffolds_refine:
      for cont in scf.contigs:
        if str(cont) not in main_strand:
            refine_contigs.add(str(cont))

    fout.write("digraph {\n")
    for edge_id, edge in edges.iteritems():
        if edge_id in main_strand:
            fout.write("{0} -> {1} [label=\"{2}\", color=\"red\"];\n".format(edge.start, edge.end, edge_id))
        elif edge_id[1:] in all_contigs:
            fout.write("{0} -> {1} [label=\"{2}\", color=\"blue\"];\n".format(edge.start, edge.end, edge_id))
        elif edge_id in refine_contigs:
            fout.write("{0} -> {1} [label=\"{2}\", color=\"green\"];\n".format(edge.start, edge.end, edge_id))
        else:
            fout.write("{0} -> {1} [label=\"{2}\"];\n".format(edge.start, edge.end, edge_id))
    fout.write("}")

def save_distance_overlap_graph(graph_file, scaffolds_in, output_file):
    max_path_len = 2 * config.ASSEMBLY_MAX_PATH_LEN
    graph, edges = ar._load_dot(graph_file)

    ordered_contigs = set()
    fout = open(output_file, "w")
    fout.write("digraph {\n")
    for scf in scaffolds_in:
        for cont in scf.contigs:
            fout.write("{0} -> {1} [label=\"{2}\", color=\"red\"];\n".format(edges[str(cont)].start, edges[str(cont)].end, str(cont)))
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
                    if found_edge in ordered_contigs and found_edge != str(prev_cont) and str(new_cont) != found_edge:
                        is_good = False
                        break

                if is_good and len(path) <= 2 * max_path_len and len(path) > 1 and (src, dst) not in mark:
                    mark.add((src, dst))
                    fout.write("{0} -> {1} [label=\"{2}\"];\n".format(src, dst, len(path) - 1))
    fout.write("}")
    fout.close()