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

logger = logging.getLogger()

#PUBLIC:
########################################################

#does the job
def refine_scaffolds(graph_file, scaffolds):
    max_path_len = config.ASSEMBLY_MAX_PATH_LEN
    logger.info("Refining with assembly graph")
    logger.debug("Max path len = {0}".format(max_path_len))
    #new_scaffolds = _insert_from_graph(graph_file, scaffolds, max_path_len)
    new_scaffolds = _insert_from_graph_experement(graph_file, scaffolds, max_path_len) #
    return new_scaffolds


#PRIVATE:
#########################################################

#ignore heavy python-graphviz
def _load_dot(filename):
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


#finds a unique path between two nodes in graph
def _get_unique_path(graph, prev_cont, next_cont, max_path_len):
    src, dst = str(prev_cont), str(next_cont)
    if not (graph.has_node(src) and graph.has_node(dst)):
        logger.debug("contigs {0} / {1} are not in the graph"
                     .format(prev_cont, next_cont))
        return None

    if graph.has_edge(src, dst):
        logger.debug("adjacent contigs {0} -- {1}".format(prev_cont, next_cont))
        return None

    paths = list(nx.all_simple_paths(graph, src, dst, max_path_len))
    if not paths:
        logger.debug("no path between {0} -- {1}".format(prev_cont, next_cont))
        return None

    if len(paths) > 1:
        logger.debug("multiple paths {0} -- {1}".format(prev_cont, next_cont))
        return None

    path = paths[0]
    path_nodes = list(map(str, path[1:-1]))
    logger.debug("unique path {0} -- {1} of length {2}"
                 .format(prev_cont, next_cont, len(path_nodes)))

    return path_nodes

def my_all_simple_paths(graph, src, dst, ordered_contigs, max_path_len):
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



def _get_unigue_path_experiment(graph, prev_cont, next_cont, max_path_len, ordered_contigs):
    src, dst =  str(prev_cont), str(next_cont)

    if not (graph.has_node(src) and graph.has_node(dst)):
        logger.debug("contigs {0} / {1} are not in the graph"
                     .format(prev_cont, next_cont))
        return None

    if graph.has_edge(src, dst):
        logger.debug("adjacent contigs {0} -- {1}".format(prev_cont, next_cont))
        return None

    paths = list(my_all_simple_paths(graph, src, dst, ordered_contigs, max_path_len))

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

def _insert_from_graph_experement(graph_file, scaffolds_in, max_path_len):
    new_scaffolds = []
    graph = _load_dot(graph_file)

    ordered_contigs = set()
    for scf in scaffolds_in:
        ordered_contigs |= set(map(lambda s: s.name, scf.contigs))

    for scf in scaffolds_in:
        new_scaffolds.append(Scaffold(scf.name))

        for prev_cont, new_cont in zip(scf.contigs[:-1], scf.contigs[1:]):
            new_scaffolds[-1].contigs.append(prev_cont)

            #find unique path
            path_nodes = _get_unigue_path_experiment(graph, prev_cont, new_cont, max_path_len, ordered_contigs)

            if not path_nodes:
                continue

            #insert contigs along the path
            for node in path_nodes:
                new_scaffolds[-1].contigs.append(Contig.from_sting(node))

        new_scaffolds[-1].contigs.append(new_cont)
    return new_scaffolds

#inserts contigs from the assembly graph into scaffolds
def _insert_from_graph(graph_file, scaffolds_in, max_path_len):
    new_scaffolds = []
    graph = _load_dot(graph_file)

    ordered_contigs = set()
    for scf in scaffolds_in:
        ordered_contigs |= set(map(lambda s: s.name, scf.contigs))

    for scf in scaffolds_in:
        new_scaffolds.append(Scaffold(scf.name))

        for prev_cont, new_cont in zip(scf.contigs[:-1], scf.contigs[1:]):
            new_scaffolds[-1].contigs.append(prev_cont)

            #find unique path
            path_nodes = _get_unique_path(graph, prev_cont, new_cont, max_path_len)

            if not path_nodes:
                continue

            #check path consistency
            consistent = True
            for node in path_nodes:
                if node[1:] in ordered_contigs:
                    logger.debug("Path inconsistency {0} -- {1}: {2}"
                                 .format(prev_cont, new_cont, node))
                    consistent = False
                    break
            if not consistent:
                continue

            #insert contigs along the path
            for node in path_nodes:
                new_scaffolds[-1].contigs.append(Contig.from_sting(node))

        new_scaffolds[-1].contigs.append(new_cont)
    return new_scaffolds
