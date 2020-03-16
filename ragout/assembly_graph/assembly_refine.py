#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module performs refinement with the assembly (overlap) graph
"""

from __future__ import absolute_import
from __future__ import division
import networkx as nx
import re
import logging
from ragout.six.moves import zip
try:
    import ragout.six.moves.queue
except ImportError:
    import queue as Queue

from ragout.shared import config
from ragout.shared.datatypes import Contig, Scaffold

logger = logging.getLogger()

def refine_scaffolds(graph_file, scaffolds, contigs_fasta):
    """
    Does the job
    """
    max_path_len = config.vals["overlap"]["max_path_len"]
    logger.info("Refining with assembly graph")
    logger.debug("Max path len = %d", max_path_len)
    graph = _load_dot(graph_file)
    _check_overaps_number(graph, contigs_fasta)
    new_scaffolds = _insert_from_graph(graph, scaffolds,
                                       max_path_len, contigs_fasta)
    _reestimate_distances(graph, new_scaffolds, contigs_fasta)
    return new_scaffolds


def _load_dot(filename):
    """
    Loads dot file (ignore heavy python-graphviz)
    """
    graph = nx.DiGraph()
    pattern = re.compile(r'"(.+)"\s*\->\s*"(.+)"\s*\[.*=.*"(.+)".*\];')
    for line in open(filename, "r").read().splitlines():
        m = pattern.match(line)
        if not m:
            continue

        v1, v2 = m.group(1), m.group(2)
        assert not graph.has_edge(v1, v2)
        graph.add_edge(v1, v2, label=m.group(3))
    return graph


def _check_overaps_number(graph, contigs_fasta):
    rate = float(len(graph.edges)) / len(contigs_fasta)
    if rate < config.vals["min_overlap_rate"]:
        logger.warning("Too few overlaps (%d) between contigs were detected "
                       "-- refine procedure will be useless. Possible reasons:"
                       "\n\n1. Some contigs output by assembler are missing\n"
                       "2. Contigs overlap not on a constant value "
                       "(like k-mer for assemblers which use debruijn graph)\n"
                       "3. Contigs ends are trimmed/postprocessed\n",
                       len(graph.edges))


def _insert_from_graph(graph, scaffolds_in, max_path_len, contigs_fasta):
    """
    Inserts contigs from the assembly graph into scaffolds
    """
    new_scaffolds = []
    ordered_contigs = set()
    for scf in scaffolds_in:
        ordered_contigs |= set([c.name() for c in scf.contigs])
    reverse_graph = graph.reverse()

    for scf in scaffolds_in:
        new_scaffolds.append(Scaffold(scf.name))

        for prev_cont, new_cont in zip(scf.contigs[:-1], scf.contigs[1:]):
            new_scaffolds[-1].contigs.append(prev_cont)

            #find contigs to insert
            path_nodes = _get_cut_vertices(graph, reverse_graph, prev_cont,
                                           new_cont, max_path_len,
                                           ordered_contigs)

            if not path_nodes:
                continue

            #insert contigs along the path
            supp_genomes = prev_cont.link.supporting_genomes
            prev_cont.link.supporting_assembly = True
            prev_cont.link.gap = config.vals["min_scaffold_gap"]
            for node in path_nodes:
                sign = 1 if node[0] == "+" else -1
                name = node[1:]

                new_contig = Contig.with_sequence(name, 
                                    len(contigs_fasta[name]), sign)
                new_contig.link.supporting_assembly = True
                new_contig.link.gap = config.vals["min_scaffold_gap"]
                new_contig.link.supporting_genomes = supp_genomes
                new_scaffolds[-1].contigs.append(new_contig)

        new_scaffolds[-1].contigs.append(scf.contigs[-1])

    return new_scaffolds


def _get_cut_vertices(graph, reverse_graph, prev_cont, next_cont,
                      max_path_len, ordered_contigs):
    """
    Finds cut vertices on a subgraph of all possible paths from one
    node to another. Corresponding contigs will be inserted into scaffolds
    between src and dst. This is a generalized version of what we have in paper
    """
    src, dst = prev_cont.signed_name(), next_cont.signed_name()

    if not (graph.has_node(src) and graph.has_node(dst)):
        #logger.debug("contigs {0} / {1} are not in the graph"
        #             .format(prev_cont, next_cont))
        return None

    if graph.has_edge(src, dst):
        #logger.debug("adjacent contigs {0} -- {1}".format(prev_cont, next_cont))
        return None

    restricted_nodes = set()
    for contig in ordered_contigs:
        restricted_nodes.add("+" + contig)
        restricted_nodes.add("-" + contig)

    induced_subgraph = _get_induced_subgraph(graph, reverse_graph, src, dst,
                                             max_path_len, restricted_nodes)

    if (not induced_subgraph.has_node(src) or
        not induced_subgraph.has_node(dst) or
        not nx.has_path(induced_subgraph, src, dst)):
        return []

    path = _shortest_path(induced_subgraph, src, dst, restricted_nodes)
    assert path is not None

    cut_vertices = set()
    for node in path[1:-1]:
        restricted_nodes.add(node)
        if (not _test_connectivity(induced_subgraph, src, dst,
                                   max_path_len, restricted_nodes)):
            cut_vertices.add(node)
        restricted_nodes.remove(node)

    ordered_cut_vertices = [p for p in path if p in cut_vertices]

    #if len(ordered_cut_vertices):
        #logger.debug("found {0} cut vertixes between {1} -- {2}"
        #             .format(len(ordered_cut_vertices), prev_cont, next_cont))

    return ordered_cut_vertices


def _get_induced_subgraph(input_graph, reverse_graph, src, dst,
                          max_path_len, restricted_nodes):
    """
    Finds subgraphs in which all possible paths between two nodes lie
    """
    def dfs(graph, vertex, end_vertex, depth, visited):
        visited.add(vertex)
        if depth == max_path_len:
            return

        for _, u in graph.edges(vertex):
            if u == end_vertex:
                visited.add(u)
                continue

            if u not in visited and u not in restricted_nodes:
                dfs(graph, u, end_vertex, depth + 1, visited)

    visited_fwd = set()
    dfs(input_graph, src, dst, 0, visited_fwd)

    visited_back = set()
    dfs(reverse_graph, dst, src, 0, visited_back)

    result = list(visited_fwd.intersection(visited_back))

    induced_digraph = nx.DiGraph()
    for node in result:
        for u, v in input_graph.edges(node):
            if v in result:
                induced_digraph.add_edge(u, v)
    return induced_digraph


def _reestimate_distances(graph, scaffolds, contigs_fasta):
    """
    Estimates distances between contigs using overlap graph
    """
    restricted_nodes = set()
    for scf in scaffolds:
        for contig in scf.contigs:
            restricted_nodes.add("+" + contig.name())
            restricted_nodes.add("-" + contig.name())

    for scf in scaffolds:
        for prev_cont, next_cont in zip(scf.contigs[:-1], scf.contigs[1:]):
            src, dst =  prev_cont.signed_name(), next_cont.signed_name()
            if graph.has_edge(src, dst):
                overlap = graph[src][dst]["label"]
                prev_cont.link.gap = -int(overlap)

            else:
                path = _shortest_path(graph, src, dst, restricted_nodes)
                if not path:
                    continue

                path_len = 0
                for node in path[1:-1]:
                    path_len += len(contigs_fasta[node[1:]])
                for n1, n2 in zip(path[:-1], path[1:]):
                    overlap = graph[n1][n2]["label"]
                    path_len -= int(overlap)

                prev_cont.link.gap = path_len


def _shortest_path(graph, src, dst, restricted_nodes):
    """
    Finds shortest path wrt to restricted nodes
    """
    queue = ragout.six.moves.queue.Queue()
    queue.put(src)
    visited = set([src])
    parent = {src : src}
    found = False
    if not graph.has_node(src):
        return None

    while not queue.empty():
        node = queue.get()

        for u in sorted(graph.neighbors(node)):
            if u == dst:
                parent[u] = node
                found = True
                break

            if u not in visited and u not in restricted_nodes:
                visited.add(u)
                queue.put(u)
                parent[u] = node

    if not found:
        return None

    path = [dst]
    cur_node = dst
    while cur_node != src:
        path.append(parent[cur_node])
        cur_node = parent[cur_node]
    return path[::-1]


def _test_connectivity(graph, start, end, max_path_len, restricted_nodes):
    """
    Quickly tests if there is a path between two nodes
    """
    class ExitSuccess(Exception):
        pass

    def dfs(node, depth):
        visited.add(node)
        if depth == max_path_len:
            return

        for _, u in graph.edges(node):
            if u == end:
                raise ExitSuccess

            if u not in visited and u not in restricted_nodes:
                dfs(u, depth + 1)

    visited = set()
    try:
        dfs(start, 0)
    except ExitSuccess:
        return True
    return False
