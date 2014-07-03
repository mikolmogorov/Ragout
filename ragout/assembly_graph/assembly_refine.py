#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module performs refinement with the assembly (overlap) graph
"""

import networkx as nx
import re
import logging
from collections import namedtuple

from ragout.shared import config
from ragout.shared.datatypes import Contig, Scaffold

logger = logging.getLogger()

def refine_scaffolds(graph_file, scaffolds, contigs_fasta):
    """
    Does the job
    """
    max_path_len = config.vals["overlap"]["max_path_len"]
    logger.info("Refining with assembly graph")
    logger.debug("Max path len = {0}".format(max_path_len))
    graph = _load_dot(graph_file)
    _check_overaps_number(graph, contigs_fasta)
    new_scaffolds = _insert_from_graph(graph, scaffolds, max_path_len)
    _reestimate_distances(graph, new_scaffolds, max_path_len, contigs_fasta)
    return new_scaffolds


def _load_dot(filename):
    """
    Loads dot file (ignore heavy python-graphviz)
    """
    graph = nx.DiGraph()
    pattern = re.compile("\"(.+)\"\s*\->\s*\"(.+)\"\s*\[.*=.*\"(.+)\".*\];")
    for line in open(filename, "r").read().splitlines():
        m = pattern.match(line)
        if not m:
            continue

        v1, v2 = m.group(1), m.group(2)
        assert not graph.has_edge(v1, v2)
        graph.add_edge(v1, v2, label=m.group(3))
    return graph


def _check_overaps_number(graph, contigs_fasta):
    rate = float(len(graph.edges())) / len(contigs_fasta)
    if rate < config.vals["min_overlap_rate"]:
        logger.warning("Too few overlaps ({0}) between contigs were detected "
                       "-- refine procedure will be useless. Possible reasons:"
                       "\n\n1. Some contigs output by assembler are missing\n"
                       "2. Contigs overlap not on a constant value "
                       "(like k-mer for assemblers which use debruijn graph)\n"
                       "3. Contigs ends are trimmed/postprocessed\n"
                       .format(len(graph.edges())))


def _insert_from_graph(graph, scaffolds_in, max_path_len):
    """
    Inserts contigs from the assembly graph into scaffolds
    """
    new_scaffolds = []
    ordered_contigs = set()
    for scf in scaffolds_in:
        ordered_contigs |= set(map(lambda s: s.name, scf.contigs))

    for scf in scaffolds_in:
        new_scaffolds.append(Scaffold(scf.name))

        for prev_cont, new_cont in zip(scf.contigs[:-1], scf.contigs[1:]):
            new_scaffolds[-1].contigs.append(prev_cont)

            #find contigs to insert
            path_nodes = _get_cut_vertices(graph, prev_cont, new_cont,
                                           max_path_len, ordered_contigs)

            if not path_nodes:
                continue

            #insert contigs along the path
            for node in path_nodes:
                sign = 1 if node[0] == "+" else -1
                name = node[1:]
                new_scaffolds[-1].contigs.append(Contig(name, sign))

        new_scaffolds[-1].contigs.append(new_cont)
    return new_scaffolds


def _reestimate_distances(graph, scaffolds, max_path_len, contigs_fasta):
    """
    Estimates distances between contigs using overlap graph
    """
    ordered_contigs = set()
    for scf in scaffolds:
        ordered_contigs |= set(map(lambda s: s.name, scf.contigs))

    for scf in scaffolds:
        for prev_cont, next_cont in zip(scf.contigs[:-1], scf.contigs[1:]):
            src, dst =  str(prev_cont), str(next_cont)
            if graph.has_edge(src, dst):
                overlap = graph[src][dst]["label"]
                prev_cont.gap = -int(overlap)

            else:
                paths = _all_simple_paths(graph, src, dst,
                                          ordered_contigs, max_path_len)
                if not paths:
                    continue
                else:
                    paths_lens = []
                    for p in paths:
                        path_len = 0
                        for node in p[1:-1]:
                            path_len += len(contigs_fasta[node[1:]])
                        for n1, n2 in zip(p[:-1], p[1:]):
                            overlap = graph[n1][n2]["label"]
                            path_len -= int(overlap)
                        paths_lens.append(path_len)
                    prev_cont.gap = _median(paths_lens)


def _get_cut_vertices(graph, prev_cont, next_cont, max_path_len,
                      ordered_contigs):
    """
    Finds cut vertices on a subgraph of all possible paths from one
    node to another. Corresponding contigs will be inserted into scaffolds
    between src and dst. This is a generalized version of what we have in paper
    """
    src, dst =  str(prev_cont), str(next_cont)

    if not (graph.has_node(src) and graph.has_node(dst)):
        logger.debug("contigs {0} / {1} are not in the graph"
                     .format(prev_cont, next_cont))
        return None

    if graph.has_edge(src, dst):
        logger.debug("adjacent contigs {0} -- {1}".format(prev_cont, next_cont))
        return None

    paths = _all_simple_paths(graph, src, dst, ordered_contigs, max_path_len)
    if not paths:
        logger.debug("no path between {0} -- {1}".format(prev_cont, next_cont))
        return None

    cut_vertices = None
    for path in paths:
        p_nodes = list(map(str, path[1:-1]))

        if not cut_vertices:
            cut_vertices = p_nodes
        else:
            cut_vertices = [p for p in cut_vertices if p in p_nodes]

    if len(cut_vertices):
        logger.debug("found {0} cut vertixes between {1} -- {2}"
                     .format(len(cut_vertices), prev_cont, next_cont))

    return cut_vertices


def _all_simple_paths(graph, src, dst, ordered_contigs, max_path_len):
    """
    Finds all possible simple paths between two nodes, that do not
    pass through the contigs (nodes) that already belong to scaffolds
    """
    answer = []

    def dfs(vertex, visit):
        for _, u in graph.edges(vertex):
            if u in visited:
                continue
            if u == dst:
                visit.append(dst)
                answer.append(list(visit))
                visit.pop()
                return

        if len(visit) >= max_path_len:
            return

        for _, u in graph.edges(vertex):
            if u not in visit and str(u)[1:] not in ordered_contigs:
                visit.append(u)
                dfs(u, visit)
                visit.pop()

    visited = [src]
    dfs(src, visited)
    return answer


def _median(values):
    """
    Not a true median, but we keep real distances
    """
    sorted_values = sorted(values)
    return sorted_values[(len(values) - 1) / 2]
