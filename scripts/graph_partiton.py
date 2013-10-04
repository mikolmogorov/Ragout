from collections import namedtuple
from breakpoint_graph import *
from itertools import product
import phylogeny as phylo
import networkx as nx
import os


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
            g.add_edge(node_1, node_2)
            g.add_edge(node_2, node_1)
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
            for neighbor in graph.neighbors(node):
                adjacencies["target"] = neighbor
                likelihood, tree = self.phylogeny.estimate_tree(adjacencies)

                if likelihood > max_likelihood:
                    max_likelihood = likelihood
                    max_vertex = neighbor
                    max_tree = tree

            g.add_edge(node, max_vertex)
            if self.debug_dir:
                debug_file = open(debug_pref + "node_{0}.dot".format(node), "w")
                phylo.tree_to_dot(max_tree, debug_file)

        return g


    def find_adjacencies(self):
        adjacencies = {}
        component_counter = 0
        for subgraph in self.graph.get_subgraphs():
            debug_pref = os.path.join(self.debug_dir, "comp{0}-".format(component_counter))
            prelinks = self.build_prelinks(subgraph, debug_pref)

            if self.debug_dir and len(subgraph) > 2:
                bg_out = open(os.path.join(self.debug_dir, "comp{0}-bg.dot".format(component_counter)), "w")
                prelinks_out = os.path.join(self.debug_dir, "comp{0}-prelinks.dot".format(component_counter))
                write_colored_dot(subgraph, bg_out)
                nx.draw_graphviz(prelinks)
                nx.write_dot(prelinks, prelinks_out)

            chosen_edges = self.split_graph(prelinks)
            for edge in chosen_edges:
                adjacencies[-edge[0]] = Connection(-edge[0], edge[1], edge[2])
                adjacencies[-edge[1]] = Connection(-edge[1], edge[0], edge[2])

            component_counter += 1

        return adjacencies


    def split_graph(self, digraph):
        max_set = None
        max_score = 0

        simple_graph = digraph.to_undirected()
        for bitset in product([0, 1], repeat=len(simple_graph.edges())):
            if not validate_set(bitset, simple_graph):
                continue

            score = 0
            for i, edge in enumerate(simple_graph.edges()):
                if bitset[i]:
                    if digraph.has_edge(edge[0], edge[1]):
                        score += 1
                    if digraph.has_edge(edge[1], edge[0]):
                        score += 1
            if score >= max_score:
                max_score = score
                max_set = bitset

        #assert max_set

        chosen_edges = [e for i, e in enumerate(simple_graph.edges()) if max_set[i]]
        ret_edges = []
        for e in chosen_edges:
            distance = self.graph.vertex_distance(e[0], e[1])
            ret_edges.append((e[0], e[1], distance))
        return ret_edges


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

