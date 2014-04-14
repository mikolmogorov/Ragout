#This module implements a breakpoint graph
#as well as the main algorithm that recovers missing 
#adjacencies
################################################

from collections import namedtuple
from itertools import chain
import os
import logging

import networkx as nx

from ragout.shared.debug import DebugConfig

Connection = namedtuple("Connection", ["start", "end"])
logger = logging.getLogger()
debugger = DebugConfig.get_instance()

#PUBLIC:
################################################


class BreakpointGraph:
    def __init__(self):
        self.bp_graph = nx.MultiGraph()
        self.targets = []
        self.references = []
        self.known_adjacencies = {}

    #builds breakpoint graph from permutations
    def build_from(self, perm_container, circular_refs):
        logger.info("Building breakpoint graph")

        for perm in perm_container.ref_perms_filtered:
            if perm.ref_id not in self.references:
                self.references.append(perm.ref_id)

        for perm in perm_container.target_perms_filtered:
            if perm.ref_id not in self.targets:
                self.targets.append(perm.ref_id)

        for perm in chain(perm_container.ref_perms_filtered,
                          perm_container.target_perms_filtered):
            circular = circular_refs if perm.ref_id in self.references else False

            prev_block = None
            for block in perm.iter_blocks(circular):
                if not prev_block:
                    prev_block = block
                    continue

                self.bp_graph.add_node(-prev_block)
                self.bp_graph.add_node(block)
                self.bp_graph.add_edge(-prev_block, block, genome_id=perm.ref_id)
                prev_block = block

    #infers missing adjacencies (the main Ragout part)
    def find_adjacencies(self, phylogeny):
        logger.info("Resolving breakpoint graph")
        chosen_edges = []
        subgraphs = nx.connected_component_subgraphs(self.bp_graph)

        for comp_id, subgraph in enumerate(subgraphs):
            trimmed_graph = self.trim_known_edges(subgraph)

            if len(trimmed_graph) < 2:
                continue

            if len(trimmed_graph) == 2:
                node_1, node_2 = trimmed_graph.nodes()
                chosen_edges.append((node_1, node_2))
                continue

            weighted_graph = self.make_weighted(trimmed_graph, phylogeny)
            matching_edges = _split_graph(weighted_graph)
            chosen_edges.extend(matching_edges)

        adjacencies = {}
        for edge in chosen_edges:
            adjacencies[-edge[0]] = Connection(-edge[0], edge[1])
            adjacencies[-edge[1]] = Connection(-edge[1], edge[0])

        if debugger.debugging:
            phylo_out = os.path.join(debugger.debug_dir, "phylogeny.txt")
            graph_out = os.path.join(debugger.debug_dir, "breakpoint_graph.dot")
            edges_out = os.path.join(debugger.debug_dir, "predicted_edges.dot")
            _output_graph(self.bp_graph, graph_out)
            _output_edges(chosen_edges, edges_out)
            _output_phylogeny(phylogeny.tree_string, self.targets[0], phylo_out)

        return adjacencies

    #removes edges with known target's adjacencies
    def trim_known_edges(self, graph):
        trimmed_graph = graph.copy()
        for v1, v2, data in graph.edges_iter(data=True):
            if not trimmed_graph.has_node(v1) or not trimmed_graph.has_node(v2):
                continue

            genome_ids = list(map(lambda e: e["genome_id"],
                                  graph[v1][v2].values()))
            target_id = self.targets[0]
            if target_id in genome_ids:
                trimmed_graph.remove_node(v1)
                trimmed_graph.remove_node(v2)

        return trimmed_graph

    #converts breakpoint graph into weighted graph
    def make_weighted(self, graph, phylogeny):
        assert len(graph) > 2
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


#PRIVATE:
###########################################################################


def _split_graph(graph):
    for v1, v2 in graph.edges_iter():
        graph[v1][v2]["weight"] = -graph[v1][v2]["weight"] #want minimum weight

    edges = nx.max_weight_matching(graph, maxcardinality=True)
    unique_edges = []
    for v1, v2 in edges.items():
        if not (v2, v1) in unique_edges:
            unique_edges.append((v1, v2))

    return unique_edges


def _update_edge(graph, v1, v2, weight):
    if not graph.has_edge(v1, v2):
        graph.add_edge(v1, v2, weight=weight)
    else:
        graph[v1][v2]["weight"] += weight

################################
#output generators

def _output_graph(graph, out_file):
    agraph = nx.write_dot(graph, out_file)


def _output_edges(edges, out_file):
    fout = open(out_file, "w")
    fout.write("graph {\n")
    for (v1, v2) in edges:
        fout.write("{0} -- {1};\n".format(v1, v2))
    fout.write("}")


def _output_phylogeny(tree_string, target_name, out_file):
    fout = open(out_file, "w")
    fout.write(tree_string + "\n")
    fout.write(target_name)
