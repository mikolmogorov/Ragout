#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module infers missing adjacencies
by recovering perfect matching
"""

from collections import namedtuple
import logging
import os

import networkx as nx

from ragout.shared.debug import DebugConfig

logger = logging.getLogger()
debugger = DebugConfig.get_instance()
Adjacency = namedtuple("Adjacency", ["block", "distance", "supporting_genomes"])


class AdjacencyInferer(object):
    def __init__(self, breakpoint_graph, phylogeny):
        self.main_graph = breakpoint_graph
        self.phylogeny = phylogeny

    def infer_adjacencies(self):
        """
        Infers missing adjacencies by recovering perfect matching
        """
        logger.info("Infering missing adjacencies")

        subgraphs = self.main_graph.connected_components()
        logger.debug("Found {0} connected components"
                             .format(len(subgraphs)))

        chosen_edges = []
        self.orphans_count = 0
        self.guessed_count = 0
        self.trimmed_count = 0
        for subgraph in subgraphs:
            chosen_edges.extend(self._process_component(subgraph))

        logger.debug("Inferred {0} adjacencies".format(len(chosen_edges)))
        logger.debug("{0} orphaned nodes".format(self.orphans_count))
        logger.debug("{0} guessed edges".format(self.guessed_count))
        logger.debug("{0} trimmed edges".format(self.trimmed_count))

        adjacencies = {}
        for node_1, node_2 in chosen_edges:
            #infinity edges correspond to joined chromosome ends -- ignore them
            if self.main_graph.is_infinity(node_1, node_2):
                continue

            distance = self.main_graph.get_distance(node_1, node_2)
            supporting_genomes = self.main_graph \
                                        .supporting_genomes(node_1, node_2)

            assert abs(node_1) != abs(node_2)
            adjacencies[node_1] = Adjacency(node_2, distance,
                                            supporting_genomes)
            adjacencies[node_2] = Adjacency(node_1, distance,
                                            supporting_genomes)

        self.main_graph.debug_output()
        self._debug_output(chosen_edges)

        return adjacencies

    def _process_component(self, subgraph):
        """
        Processes a connected component of the breakpoint graph
        """
        adjacency = subgraph.to_weighted_graph(self.phylogeny)
        trimmed_graph = self._trim_known_edges(adjacency)
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

            matching_edges = _min_weight_matching(trim_subgraph)

            for edge in matching_edges:
                for n in edge:
                    unused_nodes.remove(n)

            chosen_edges.extend(matching_edges)

        #predicting target-specific rearrangement
        #NO!
        #if len(unused_nodes) == 2:
        #    node_1, node_2 = tuple(unused_nodes)
        #    cycle = subgraph.alternating_cycle(node_1, node_2)
        #    if (abs(node_1) != abs(node_2) and cycle in [2, 3]):
        #        self.guessed_count += 1
        #        chosen_edges.append((node_1, node_2))
        #        unused_nodes.clear()
        self.orphans_count += len(unused_nodes)

        return chosen_edges

    def _trim_known_edges(self, graph):
        """
        Removes edges with known target adjacencies (red edges from paper)
        """
        trimmed_graph = graph.copy()
        for v1, v2 in graph.edges_iter():
            if not trimmed_graph.has_node(v1) or not trimmed_graph.has_node(v2):
                continue

            genome_ids = self.main_graph.supporting_genomes(v1, v2)
            if self.main_graph.target in genome_ids:
                for node in [v1, v2]:
                    trimmed_graph.remove_node(node)
                self.trimmed_count += 1

        return trimmed_graph

    def _debug_output(self, chosen_edges):
        if not debugger.debugging:
            return

        phylo_out = os.path.join(debugger.debug_dir, "phylogeny.txt")
        edges_out = os.path.join(debugger.debug_dir, "predicted_edges.dot")
        _output_edges(chosen_edges, edges_out)
        _output_phylogeny(self.phylogeny.tree_string, self.main_graph.target,
                          phylo_out)


def _min_weight_matching(graph):
    """
    Finds a perfect matching with minimum weight
    """
    for v1, v2 in graph.edges_iter():
        graph[v1][v2]["weight"] = -graph[v1][v2]["weight"] #want minimum weght

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
