from collections import namedtuple
from breakpoint_graph import *
from itertools import product
from networkx import Graph


Connection = namedtuple("Connection", ["start", "end", "distance"])


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


def split_graph(graph):
    max_set = None
    max_score = 0
    max_edges = 0
    for bitset in product([0, 1], repeat=len(graph.edges())):
        if not validate_set(bitset, graph):
            continue

        score = 0
        num_edges = 0
        for i, edge in enumerate(graph.edges()):
            if bitset[i]:
                num_edges += 1
                score += 1
        if num_edges > max_edges:
            max_edges, max_score = num_edges, score
            max_set = bitset
        elif num_edges == max_edges and score > max_score:
            max_score = score
            max_set = bitset

    assert max_set
    filtered = [e for i, e in enumerate(graph.edges()) if max_set[i]]
    ret_edges = [(e[0], e[1], graph[e[0]][e[1]]["distance"]) for e in filtered]
    return ret_edges


class AdjacencyFinder:
    def __init__(self, graph, sibelia_output):
        self.graph = graph
        self.sibelia_output = sibelia_output


    def compress_graph(self, graph):
        g = Graph()
        for u, v, data in graph.edges_iter(data=True):
            weight = -data["color"]
            if g.has_edge(u, v):
                g[u][v]["weight"] = max(weight, g[u][v]["weight"])
            else:
                distance = self.graph.vertex_distance(u, v)
                g.add_edge(u, v, weight=weight, distance=distance)
        return g


    def find_adjacencies(self):
        adjacencies = {}
        for subgraph in self.graph.get_unresolved_subgraphs():
            compressed = self.compress_graph(subgraph)
            edges = split_graph(compressed)
            for edge in edges:
                adjacencies[-edge[0]] = Connection(-edge[0], edge[1], edge[2])
                adjacencies[-edge[1]] = Connection(-edge[1], edge[0], edge[2])

        return adjacencies
