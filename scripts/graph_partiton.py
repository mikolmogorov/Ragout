from collections import namedtuple
from breakpoint_graph import *
from itertools import product
import phylogeny as phylo
import networkx as nx
import os
import math


Connection = namedtuple("Connection", ["start", "end", "distance"])


class AdjacencyFinder:
    def __init__(self, graph, phylogeny, debug_dir=None):
        self.graph = graph
        self.phylogeny = phylogeny
        self.debug_dir = debug_dir


    def build_prelinks(self, graph, debug_pref):
        g = nx.DiGraph()
        g.add_nodes_from(graph.nodes())

        #trivial case
        if len(graph) == 2:
            node_1, node_2 = graph.nodes()
            g.add_edge(node_1, node_2, weight=1)
            g.add_edge(node_2, node_1, weight=1)
            return g

        #non-trivial
        for node in graph.nodes():
            adjacencies = {}
            for neighbor in graph.neighbors(node):
                for edge in graph[node][neighbor].values():
                    adjacencies[edge["ref_id"]] = neighbor

            max_likelihood = float("-inf")
            max_vertex = None
            max_tree = None
            min_breaks = 99999

            edges_to_add = []

            for neighbor in graph.neighbors(node):
                adjacencies["target"] = neighbor
                likelihood, nbreaks, tree = self.phylogeny.estimate_tree(adjacencies)

                if likelihood > max_likelihood:
                    max_likelihood = likelihood
                    max_vertex = neighbor
                    max_tree = tree

                #if nbreaks < min_breaks:
                #    edges_to_add = [(neighbor, likelihood)]
                #    min_breaks = nbreaks
                #elif nbreaks == min_breaks:
                #    edges_to_add.append((neighbor, likelihood))

                g.add_edge(node, neighbor, weight=likelihood)

            #for edge in edges_to_add:
            #    g.add_edge(node, edge[0], weight=edge[1])

            #g.add_edge(node, max_vertex)
            if self.debug_dir:
                debug_file = open(debug_pref + "node_{0}.dot".format(node), "w")
                phylo.tree_to_dot(max_tree, debug_file)

        return g


    def find_adjacencies(self):
        adjacencies = {}
        component_counter = 0
        for subgraph in self.graph.get_subgraphs():
            debug_pref = os.path.join(self.debug_dir, "comp{0}-".format(component_counter))
            digraph = self.build_prelinks(subgraph, debug_pref)

            weighted_graph = nx.Graph()
            min_weight = float("+inf")
            for edge in digraph.edges_iter(data=True):
                weight = edge[2]["weight"]
                if not weighted_graph.has_edge(edge[0], edge[1]):
                    weighted_graph.add_edge(edge[0], edge[1], weight=weight)
                else:
                    weighted_graph[edge[0]][edge[1]]["weight"] += weight
                    min_weight = min(weighted_graph[edge[0]][edge[1]]["weight"], min_weight)

            for edge in weighted_graph.edges_iter():
                weighted_graph[edge[0]][edge[1]]["weight"] -= min_weight
                weighted_graph[edge[0]][edge[1]]["label"] = "{0:5.2f}".format(weighted_graph[edge[0]][edge[1]]["weight"])
                #print weighted_graph[edge[0]][edge[1]]["label"]

            chosen_edges = self.split_graph(weighted_graph)
            #edges = nx.max_weight_matching(weighted_graph)
            #print chosen_edges
            #print edges

            if self.debug_dir and len(subgraph) > 2:
                bg_out = open(os.path.join(self.debug_dir, "comp{0}-bg.dot".format(component_counter)), "w")
                prelinks_out = os.path.join(self.debug_dir, "comp{0}-prelinks.dot".format(component_counter))
                write_colored_dot(subgraph, bg_out)
                nx.draw_graphviz(weighted_graph)
                nx.write_dot(weighted_graph, prelinks_out)

            for edge in chosen_edges:
                adjacencies[-edge[0]] = Connection(-edge[0], edge[1], edge[2])
                adjacencies[-edge[1]] = Connection(-edge[1], edge[0], edge[2])

            component_counter += 1

        return adjacencies


    def split_graph(self, graph):
        edges = nx.max_weight_matching(graph, maxcardinality=True)
        unique_edges = []
        for v1, v2 in edges.iteritems():
            if not (v2, v1) in unique_edges:
                unique_edges.append((v1, v2))

        print unique_edges
        ret_edges = []
        for e in unique_edges:
            distance = self.graph.vertex_distance(e[0], e[1])
            ret_edges.append((e[0], e[1], distance))
        return ret_edges


    """
    def split_graph(self, graph):
        max_set = None
        max_score = 0

        counter = 0
        simple_graph = graph
        #simple_graph = digraph.to_undirected()
        #simple_graph = nx.Graph()
        #for edge in digraph.edges_iter(data=True):
        #    weight = edge[2]["weight"]
        #    if not simple_graph.has_edge(edge[0], edge[1]):
        #        simple_graph.add_edge(edge[0], edge[1], label=weight)
        #    else:
        #        simple_graph[edge[0]][edge[1]]["label"] += weight

        #if self.debug_dir and len(subgraph) > 2:
            #prelinks_out = os.path.join(self.debug_dir, "comp{0}-prelinks.dot".format(component_counter))
            #nx.write_dot(prelinks, prelinks_out)

        for bitset in product([0, 1], repeat=len(simple_graph.edges())):
            counter += 1
            #print counter, 2 ** len(simple_graph.edges())
            if not validate_set(bitset, simple_graph):
                continue

            score = 0.0
            #score = 0
            for i, edge in enumerate(simple_graph.edges()):
                if bitset[i]:
                    #if digraph.has_edge(edge[0], edge[1]):
                        #score += digraph[edge[0]][edge[1]]["weight"]
                        #score += 1
                    #if digraph.has_edge(edge[1], edge[0]):
                        #score += 1
                        #score += digraph[edge[1]][edge[0]]["weight"]
                    score += simple_graph[edge[0]][edge[1]]["label"]
            if score >= max_score:
                max_score = score
                max_set = bitset

        #assert max_set

        chosen_edges = [e for i, e in enumerate(simple_graph.edges()) if max_set[i]]
        print chosen_edges
        ret_edges = []
        for e in chosen_edges:
            distance = self.graph.vertex_distance(e[0], e[1])
            ret_edges.append((e[0], e[1], distance))
        return ret_edges
"""


def validate_set(bitset, graph):
    used_vertex = set()
    for i, edge in enumerate(graph.edges()):
        if not bitset[i]:
            continue

        if edge[0] in used_vertex or edge[1] in used_vertex:
            return False
        else:
            used_vertex.add(edge[0])
            used_vertex.add(edge[1])
    return True

