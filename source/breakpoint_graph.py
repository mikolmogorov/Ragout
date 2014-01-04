import networkx as nx
from collections import namedtuple
import os
import logging

from permutation import *
from debug import DebugConfig, write_dot
import phylogeny as phylo


Connection = namedtuple("Connection", ["start", "end"])
logger = logging.getLogger()

class BreakpointGraph:
    def __init__(self):
        self.bp_graph = nx.MultiGraph()
        self.targets = []
        self.references = []
        self.known_adjacencies = {}


    def build_from(self, perm_container, circular):
        logger.info("Building breakpoint graph")

        for perm in perm_container.ref_perms_filtered:
            if perm.ref_id not in self.references:
                self.references.append(perm.ref_id)

            prev = None
            for block in perm.iter_blocks(circular):
                if not prev:
                    prev = block
                    continue

                left_block = prev
                right_block = block
                self.bp_graph.add_node(-left_block)
                self.bp_graph.add_node(right_block)
                self.bp_graph.add_edge(-left_block, right_block, ref_id=perm.ref_id)

                prev = block

        for perm in perm_container.target_perms_filtered:
            if perm.ref_id not in self.targets:
                self.targets.append(perm.ref_id)

            prev = None
            for block in perm.iter_blocks(False):
                if not prev:
                    prev = block
                    continue

                self.known_adjacencies[-prev] = block
                self.known_adjacencies[block] = -prev
                prev = block


    def find_adjacencies(self, phylogeny):
        logger.info("Resolving breakpoint graph")
        chosen_edges = []
        subgraphs = nx.connected_component_subgraphs(self.bp_graph)


        for comp_id, subgraph in enumerate(subgraphs):
            known_adjacencies, trimmed_graph = self.trim_known_edges(subgraph)
            chosen_edges.extend(known_adjacencies)

            if len(trimmed_graph) < 2:
                continue

            if len(trimmed_graph) == 2:
                node_1, node_2 = trimmed_graph.nodes()
                chosen_edges.append((node_1, node_2))
                continue

            weighted_graph = self.make_weighted(trimmed_graph, phylogeny)
            matching_edges = split_graph(weighted_graph)
            chosen_edges.extend(matching_edges)

            if DebugConfig.get_writer().debugging:
                debug_dir = DebugConfig.get_writer().debug_dir
                debug_draw_component(comp_id, weighted_graph, subgraph, debug_dir)

        adjacencies = {}
        for edge in chosen_edges:
            adjacencies[-edge[0]] = Connection(-edge[0], edge[1])
            adjacencies[-edge[1]] = Connection(-edge[1], edge[0])

        return adjacencies


    def trim_known_edges(self, graph):
        known_edges = []
        trimmed_graph = graph.copy()
        for v1, v2, data in graph.edges_iter(data=True):
            if not trimmed_graph.has_node(v1) or not trimmed_graph.has_node(v2):
                continue

            if self.known_adjacencies.get(v1, None) == v2:
                trimmed_graph.remove_node(v1)
                trimmed_graph.remove_node(v2)
                known_edges.append((v1, v2))

        return known_edges, trimmed_graph


    def make_weighted(self, graph, phylogeny):
        assert len(graph) > 2
        g = nx.Graph()
        g.add_nodes_from(graph.nodes())
        target_id = self.targets[0]

        for node in graph.nodes():
            adjacencies = {}
            for neighbor in graph.neighbors(node):
                for edge in graph[node][neighbor].values():
                    adjacencies[edge["ref_id"]] = neighbor

            for ref_id in self.references:
                if ref_id not in adjacencies:
                    adjacencies[ref_id] = None  #"void" state in paper

            for neighbor in graph.neighbors(node):
                break_weight = 0.0
                if not (self.known_adjacencies.get(node, None) == neighbor):
                    adjacencies[target_id] = neighbor
                    break_weight = phylogeny.estimate_tree(adjacencies)

                update_edge(g, node, neighbor, break_weight)

        return g


def split_graph(graph):
    for v1, v2 in graph.edges_iter():
        graph[v1][v2]["weight"] = -graph[v1][v2]["weight"] #want minimum weight

    edges = nx.max_weight_matching(graph, maxcardinality=True)
    unique_edges = []
    for v1, v2 in edges.iteritems():
        if not (v2, v1) in unique_edges:
            unique_edges.append((v1, v2))

    return unique_edges


def update_edge(graph, v1, v2, weight):
    if not graph.has_edge(v1, v2):
        graph.add_edge(v1, v2, weight=weight)
    else:
        graph[v1][v2]["weight"] += weight


def debug_draw_component(comp_id, weighted_graph, breakpoint_graph, debug_dir):
    if len(breakpoint_graph) == 2:
        return

    for e in weighted_graph.edges_iter():
        weighted_graph[e[0]][e[1]]["label"] = ("{0:7.4f}"
                                    .format(weighted_graph[e[0]][e[1]]["weight"]))
    bg_out = os.path.join(debug_dir, "comp{0}-bg.dot".format(comp_id))
    weighted_out = os.path.join(debug_dir, "comp{0}-weighted.dot".format(comp_id))
    write_dot(breakpoint_graph, open(bg_out, "w"))
    write_dot(weighted_graph, open(weighted_out, "w"))
