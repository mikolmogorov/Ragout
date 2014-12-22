#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module implements a breakpoint graph
as well as the main algorithm which recovers missing
adjacencies
"""

from collections import namedtuple
from itertools import chain
import os
import logging

import networkx as nx

from ragout.shared.debug import DebugConfig

Adjacency = namedtuple("Adjacency", ["block", "distance", "supporting_genomes"])
logger = logging.getLogger()
debugger = DebugConfig.get_instance()


class BreakpointGraph:
    """
    Breakpoint graph implementation, as it is written in paper
    """
    def __init__(self):
        self.bp_graph = nx.MultiGraph()
        self.targets = []
        self.references = []
        self.known_adjacencies = {}

    def build_from(self, perm_container, recipe):
        """
        Builds breakpoint graph from permutations
        """
        logger.debug("Building breakpoint graph")

        for perm in perm_container.ref_perms:
            if perm.genome_name not in self.references:
                self.references.append(perm.genome_name)

        for perm in perm_container.target_perms:
            if perm.genome_name not in self.targets:
                self.targets.append(perm.genome_name)

        for perm in chain(perm_container.ref_perms,
                          perm_container.target_perms):

            if len(perm.blocks) < 2:
                continue

            for prev_block, next_block in perm.iter_pairs():
                self.bp_graph.add_node(-prev_block.signed_id())
                self.bp_graph.add_node(next_block.signed_id())

                distance = next_block.start - prev_block.end
                assert distance >= 0
                self.bp_graph.add_edge(-prev_block.signed_id(),
                                       next_block.signed_id(),
                                       genome_id=perm.genome_name,
                                       distance=distance)

            if (perm.genome_name in self.references and
                not recipe["genomes"][perm.genome_name]["draft"]):

                distance = (perm.chr_len - perm.blocks[-1].end +
                            perm.blocks[0].start)
                assert distance >= 0

                if recipe["genomes"][perm.genome_name]["circular"]:
                    self.bp_graph.add_edge(-perm.blocks[-1].signed_id(),
                                           perm.blocks[0].signed_id(),
                                           genome_id=perm.genome_name,
                                           distance=distance)
                else:
                    self.bp_graph.add_edge(-perm.blocks[-1].signed_id(),
                                           perm.blocks[0].signed_id(),
                                           genome_id=perm.genome_name,
                                           distance=distance,
                                           infinity=True)

        logger.debug("Built graph with {0} nodes".format(len(self.bp_graph)))

    def find_adjacencies(self, phylogeny):
        """
        Infers missing adjacencies (the main Ragout part)
        """
        logger.info("Resolving breakpoint graph")

        subgraphs = list(nx.connected_component_subgraphs(self.bp_graph))
        logger.debug("Found {0} connected components"
                     .format(len(subgraphs)))

        chosen_edges = []
        self.orphans_count = 0
        self.guessed_count = 0
        self.trimmed_count = 0
        for subgraph in subgraphs:
            chosen_edges.extend(self._process_component(subgraph, phylogeny))

        logger.debug("Inferred {0} adjacencies".format(len(chosen_edges)))
        logger.debug("{0} orphaned nodes".format(self.orphans_count))
        logger.debug("{0} guessed edges".format(self.guessed_count))
        logger.debug("{0} trimmed edges".format(self.trimmed_count))

        adjacencies = {}
        for edge in chosen_edges:
            #infinity edges correspond to joined chromosome ends -- ignore them
            if self._is_infinity(edge[0], edge[1]):
                continue

            distance = self._get_distance(edge[0], edge[1])
            supporting_genomes = []
            if self.bp_graph.has_edge(edge[0], edge[1]):
                for e in self.bp_graph[edge[0]][edge[1]].values():
                    supporting_genomes.append(e["genome_id"])
            adjacencies[-edge[0]] = Adjacency(edge[1], distance,
                                              supporting_genomes)
            adjacencies[-edge[1]] = Adjacency(edge[0], distance,
                                              supporting_genomes)

        if debugger.debugging:
            phylo_out = os.path.join(debugger.debug_dir, "phylogeny.txt")
            graph_out = os.path.join(debugger.debug_dir, "breakpoint_graph.dot")
            edges_out = os.path.join(debugger.debug_dir, "predicted_edges.dot")
            _output_graph(self.bp_graph, graph_out)
            _output_edges(chosen_edges, edges_out)
            _output_phylogeny(phylogeny.tree_string, self.targets[0], phylo_out)

        return adjacencies

    def _process_component(self, subgraph, phylogeny):
        """
        Processes a connected component of the breakpoint graph
        """
        weighted_graph = self._make_weighted(subgraph, phylogeny)
        trimmed_graph = self._trim_known_edges(weighted_graph)
        unused_nodes = set(trimmed_graph.nodes())

        chosen_edges = []
        for trim_subgraph in nx.connected_component_subgraphs(trimmed_graph):
            if len(trim_subgraph) < 2:
                continue

            if len(trim_subgraph) == 2:
                chosen_edges.append(tuple(trim_subgraph.nodes()))
                for n in trim_subgraph.nodes():
                    unused_nodes.remove(n)
                continue

            matching_edges = _split_graph(trim_subgraph)

            for edge in matching_edges:
                for n in edge:
                    unused_nodes.remove(n)
            chosen_edges.extend(matching_edges)

        #check if there are only 2 nodes left
        #if len(unused_nodes) == 2:
        #    self.guessed_count += 1
        #    chosen_edges.append(tuple(unused_nodes))
        #    unused_nodes.clear()
        self.orphans_count += len(unused_nodes)

        return chosen_edges


    def _trim_known_edges(self, graph):
        """
        Removes edges with known adjacencies in target (red edges from paper)
        """
        trimmed_graph = graph.copy()
        for v1, v2, data in graph.edges_iter(data=True):
            if not trimmed_graph.has_node(v1) or not trimmed_graph.has_node(v2):
                continue

            genome_ids = list(map(lambda e: e["genome_id"],
                                  self.bp_graph[v1][v2].values()))
            target_id = self.targets[0]
            if target_id in genome_ids:
                for node in [v1, v2]:
                    trimmed_graph.remove_node(node)
                self.trimmed_count += 1


        return trimmed_graph

    def _make_weighted(self, graph, phylogeny):
        """
        Converts a breakpoint graph into a weighted graph
        """
        assert len(graph) >= 2
        g = nx.Graph()
        g.add_nodes_from(graph.nodes())
        target_id = self.targets[0]

        for node in graph.nodes():
            adjacencies = {}
            for neighbor in graph.neighbors(node):
                for edge in graph[node][neighbor].values():
                    adjacencies[edge["genome_id"]] = neighbor

            for ref_id in self.references:
                if ref_id not in adjacencies:
                    adjacencies[ref_id] = None  #"void" state in paper

            for neighbor in graph.neighbors(node):
                adjacencies[target_id] = neighbor
                break_weight = phylogeny.estimate_tree(adjacencies)

                _update_edge(g, node, neighbor, break_weight)

        return g

    def _is_infinity(self, node_1, node_2):
        if not self.bp_graph.has_edge(node_1, node_2):
            return False

        for edge_data in self.bp_graph[node_1][node_2].values():
            if "infinity" in edge_data:
                return True
        return False

    def _get_distance(self, node_1, node_2):
        """
        Tries to guess the distance between synteny blocks
        in a target genome
        """
        DEFAULT_DISTANCE = 0
        if not self.bp_graph.has_edge(node_1, node_2):
            return DEFAULT_DISTANCE
        distances = [e["distance"]
                     for e in self.bp_graph[node_1][node_2].values()]
        return _median(distances) #currently, just a median :(


def _split_graph(graph):
    """
    Finds a perfect matching with minimum weight
    """
    for v1, v2 in graph.edges_iter():
        graph[v1][v2]["weight"] = -graph[v1][v2]["weight"] #want minimum weight

    MIN_LOG_SIZE = 20
    if len(graph) > MIN_LOG_SIZE:
        logger.debug("Finding perfect matching for a component of "
                     "size {0}".format(len(graph)))
    edges = nx.max_weight_matching(graph, maxcardinality=True)
    unique_edges = set()
    for v1, v2 in edges.items():
        if not (v2, v1) in unique_edges:
            unique_edges.add((v1, v2))

    return list(unique_edges)


def _update_edge(graph, v1, v2, weight):
    """
    Helper function to update edge's weight
    """
    if not graph.has_edge(v1, v2):
        graph.add_edge(v1, v2, weight=weight)
    else:
        graph[v1][v2]["weight"] += weight


def _median(values):
    """
    Not a true median, but we keep real distances
    """
    sorted_values = sorted(values)
    return sorted_values[(len(values) - 1) / 2]


def _output_graph(graph, out_file):
    """
    Outputs graph in dot format
    """
    with open(out_file, "w") as fout:
        fout.write("graph {\n")
        for v1, v2, data in graph.edges_iter(data=True):
            fout.write("{0} -- {1}".format(v1, v2))
            if len(data):
                extra = list(map(lambda (k, v) : "{0}=\"{1}\"".format(k, v),
                                 data.items()))
                fout.write(" [" + ", ".join(extra) + "]")
            fout.write(";\n")
        fout.write("}")


def _output_edges(edges, out_file):
    """
    Outputs list of edges in dot format
    """
    with open(out_file, "w") as fout:
        fout.write("graph {\n")
        for (v1, v2) in edges:
            fout.write("{0} -- {1};\n".format(v1, v2))
        fout.write("}")


def _output_phylogeny(tree_string, target_name, out_file):
    """
    Outputs phylogenetic tree in plain text
    """
    with open(out_file, "w") as fout:
        fout.write(tree_string + "\n")
        fout.write(target_name)
