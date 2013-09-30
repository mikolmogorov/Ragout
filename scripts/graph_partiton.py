from collections import namedtuple
from breakpoint_graph import *
from itertools import product
from networkx import Graph
from phylogeny import Phylogeny


Connection = namedtuple("Connection", ["start", "end", "distance"])


class AdjacencyFinder:
    def __init__(self, graph, phylogeny):
        self.graph = graph
        self.phylogeny = phylogeny


    def compress_graph(self, graph):
        g = Graph()

        #trivial case
        if len(graph) == 2:
            node_1, node_2 = graph.nodes()
            g.add_nodes_from(graph.nodes())
            distance = self.graph.vertex_distance(node_1, node_2)
            g.add_edge(node_1, node_2, distance=distance)
            g.node[node_1]["pair"] = node_2
            g.node[node_2]["pair"] = node_1
            return g

        #non-trivial
        for node in graph.nodes():
            adjacencies = {}
            for neighbor in graph.neighbors(node):
                distance = self.graph.vertex_distance(node, neighbor)
                g.add_edge(node, neighbor, distance=distance)
                #print graph[node][neighbor]
                for edge in graph[node][neighbor].values():
                    adjacencies[edge["ref_id"]] = neighbor

            max_likelihood = float("-inf")
            max_vertex = None
            for neighbor in graph.neighbors(node):
                adjacencies["target"] = neighbor
                likelihood = self.phylogeny.estimate_tree(adjacencies)
                if likelihood > max_likelihood:
                    max_likelihood = likelihood
                    max_vertex = neighbor

            g.node[node]["pair"] = max_vertex

        return g


    def find_adjacencies(self):
        adjacencies = {}
        for subgraph in self.graph.get_unresolved_subgraphs():
            compressed = self.compress_graph(subgraph)
            #print "==="
            edges = split_graph(compressed)
            for edge in edges:
                #print edge
                adjacencies[-edge[0]] = Connection(-edge[0], edge[1], edge[2])
                adjacencies[-edge[1]] = Connection(-edge[1], edge[0], edge[2])

        return adjacencies


def split_graph(graph):
    max_set = None
    max_score = 0
    max_edges = 0
    for bitset in product([0, 1], repeat=len(graph.edges())):
        if not validate_set(bitset, graph):
            continue

        score = 0
        for i, edge in enumerate(graph.edges()):
            if bitset[i]:
                if graph.node[edge[0]]["pair"] == edge[1]:
                    score += 1
                if graph.node[edge[1]]["pair"] == edge[0]:
                    score += 1
        if score > max_score:
            max_score = score
            max_set = bitset

    assert max_set
    filtered = [e for i, e in enumerate(graph.edges()) if max_set[i]]
    ret_edges = [(e[0], e[1], graph[e[0]][e[1]]["distance"]) for e in filtered]
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

