#This module performs refinement with the assembly graph
#########################################################

import networkx as nx
import re
import logging
from collections import namedtuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ragout.shared import config
from ragout.shared.datatypes import Contig, Scaffold

Edge = namedtuple("Edge", ["start", "end"])
logger = logging.getLogger()

#PUBLIC:
########################################################

#does the job
def refine_scaffolds(graph_file, scaffolds):
    max_path_len = config.ASSEMBLY_MAX_PATH_LEN
    logger.info("Refining with assembly graph")
    logger.debug("Max path len = {0}".format(max_path_len))
    new_scaffolds = _insert_from_graph(graph_file, scaffolds, max_path_len)
    return new_scaffolds


#PRIVATE:
#########################################################

#ignore heavy python-graphviz
def _load_dot(filename):
    graph = nx.MultiDiGraph()
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


#check if there is no multiedges along the path
def _check_unique(graph, path):
    for v1, v2 in zip(path[:-1], path[1:]):
        assert graph.has_edge(v1, v2)
        if len(graph[v1][v2]) > 1:
            return False
    return True


#finds a unique path between two nodes in graph
def _get_unique_path(graph, edges, prev_cont, new_cont, max_path_len):
    try:
        src = edges[str(prev_cont)].end
        dst = edges[str(new_cont)].start
    except KeyError:
        logger.debug("contigs are not in the graph")
        return None

    if src == dst:
        logger.debug("adjacent contigs {0} -- {1}".format(prev_cont, new_cont))
        return None

    if not nx.has_path(graph, src, dst):
        logger.debug("no path {0} -- {1}".format(prev_cont, new_cont))
        return None

    paths = [p for p in nx.all_shortest_paths(graph, src, dst)]
    if len(paths) > 1 or not _check_unique(graph, paths[0]):
        logger.debug("multiple paths {0} -- {1}".format(prev_cont, new_cont))
        return None

    path = paths[0]
    if len(path) > max_path_len:
        logger.debug("too long path {0} -- {1} of length {2}"
                       .format(prev_cont, new_cont, len(path)))
        return None

    path_edges = []
    for p_start, p_end in zip(path[:-1], path[1:]):
        found_edge = None
        for edge_id, edge in edges.items():
            if edge == Edge(p_start, p_end):
                found_edge = edge_id
                break
        assert found_edge
        path_edges.append(found_edge)

    logger.debug("unique path {0} -- {1} of length {2}"
                 .format(prev_cont, new_cont, len(path)))

    return path_edges


#inserts contigs from the assembly graph into scaffolds
def _insert_from_graph(graph_file, scaffolds_in, max_path_len):
    new_scaffolds = []
    graph, edges = _load_dot(graph_file)

    ordered_contigs = set()
    for scf in scaffolds_in:
        ordered_contigs |= set(map(lambda s: s.name, scf.contigs))

    for scf in scaffolds_in:
        new_scaffolds.append(Scaffold(scf.name))

        for prev_cont, new_cont in zip(scf.contigs[:-1], scf.contigs[1:]):
            new_scaffolds[-1].contigs.append(prev_cont)

            #find unique path
            path_edges = _get_unique_path(graph, edges, prev_cont, new_cont, max_path_len)
            if path_edges is None:
                continue

            #check path consistency
            consistent = True
            for edge in path_edges:
                if edge[1:] in ordered_contigs:
                    logger.debug("Path inconsistency {0} -- {1}: {2}"
                                 .format(prev_cont, new_cont, edge))
                    consistent = False
                    break
            if not consistent:
                continue

            #insert contigs along the path
            for edge in path_edges:
                new_scaffolds[-1].contigs[-1].gap = 0
                new_scaffolds[-1].contigs.append(Contig.from_sting(edge))
                new_scaffolds[-1].contigs[-1].gap = 0

        new_scaffolds[-1].contigs.append(new_cont)
    return new_scaffolds
