"""
This module performs refinement with the assembly graph
"""

import networkx as nx
import re
import logging
from collections import namedtuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ragout.shared import config
from ragout.shared.datatypes import Contig, Scaffold

logger = logging.getLogger()

def refine_scaffolds(graph_file, scaffolds):
    """
    Does the job
    """
    max_path_len = config.ASSEMBLY_MAX_PATH_LEN
    logger.info("Refining with assembly graph")
    logger.debug("Max path len = {0}".format(max_path_len))
    graph = _load_dot(graph_file)
    new_scaffolds = _insert_from_graph(graph, scaffolds, max_path_len)
    _reestimate_distances(graph, new_scaffolds, max_path_len)
    return new_scaffolds


def _load_dot(filename):
    """
    Loads dot file (ignore heavy python-graphviz)
    """
    graph = nx.MultiDiGraph()
    pattern = re.compile("\"(.+)\"\s*\->\s*\"(.+)\"\s*\[.*=.*\"(.+)\".*\];")
    for line in open(filename, "r").read().splitlines():
        m = pattern.match(line)
        if not m:
            continue

        v1, v2 = m.group(1), m.group(2)
        graph.add_node(v1)
        graph.add_node(v2)
        graph.add_edge(v1, v2, label=m.group(3))
    return graph


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

            #find unique path
            path_nodes = _get_unigue_path_experiment(graph, prev_cont, new_cont,
                                                     max_path_len, ordered_contigs)

            if not path_nodes:
                continue

            #insert contigs along the path
            for node in path_nodes:
                new_scaffolds[-1].contigs.append(Contig.from_sting(node))

        new_scaffolds[-1].contigs.append(new_cont)
    return new_scaffolds


def _reestimate_distances(graph, scaffolds, max_path_len):
    """
    Estimates distances between contigs using overlap graph
    """
    pass


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

    path_nodes = None
    for path in paths:
        p_nodes = list(map(str, path[1:-1]))

        if not path_nodes:
            path_nodes = p_nodes
        else:
            for node in path_nodes:
                if p_nodes.count(node) == 0:
                    path_nodes.remove(node)

    if path_nodes:
        logger.debug("unique path {0} -- {1} of length {2}"
                 .format(prev_cont, next_cont, len(path_nodes)))

    return path_nodes


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
