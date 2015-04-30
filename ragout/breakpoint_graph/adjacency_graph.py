#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
Adjacency graph -- weighted graph with black edges.
Provides some useful methods for analysis
"""

from collections import defaultdict

import networkx as nx

from ragout.shared.priority_queue import PriorityQueue

class AdjacencyGraph(object):
    """
    A wrapper for convenient use
    """
    def __init__(self, bp_graph, phylogeny):
        self.graph = self._from_bp_graph(bp_graph, phylogeny)
        self.orphan_weight = 0

    def _from_bp_graph(self, bp_graph, phylogeny):
        """
        Constructs adjacency graph (weighted with black edges)
        """
        weighted_graph = bp_graph.to_weighted_graph(phylogeny)
        adj_graph = nx.MultiGraph()
        adj_graph.add_edges_from(weighted_graph.edges(data=True))
        black_nodes = set()
        for node_1, node_2 in bp_graph.contig_ends:
            adj_graph.add_edge(node_1, node_2, black=True)
            black_nodes.add(node_1)
            black_nodes.add(node_2)

        for node in adj_graph.nodes():
            if node not in black_nodes:
                adj_graph.remove_node(node)
        return adj_graph

    def _neighbors(self, node, colored):
        neighbors = []
        for neighbor in self.graph[node]:
            for edge in self.graph[node][neighbor].values():
                if ("black" not in edge) == colored:
                    neighbors.append(neighbor)
        return neighbors

    def _edge_weight(self, node_1, node_2, colored):
        if colored:
            if self.graph.has_edge(node_1, node_2):
                for value in self.graph[node_1][node_2].values():
                    if "weight" in value:
                        return value["weight"]
            return self.orphan_weight
        else:
            return 0

    def shortest_path(self, src, dst, prohibited_nodes, orphaned_nodes):
        """
        Finds shortest path wrt to prohibited and orphaned nodes
        """
        #logger.debug("Finding path from {0} to {1}".format(src, dst))
        dist = defaultdict(lambda: float("inf"))
        dist[src] = 0
        parent = {}
        queue = PriorityQueue()
        queue.insert((src, True), 0)

        found = False
        while queue.get_length():
            cur_dist, (cur_node, colored) = queue.pop()
            if cur_node == dst:
                if not colored:
                    found = True
                    break
                else:
                    continue

            if cur_node != src and cur_node in prohibited_nodes:
                continue

            neighbors = self._neighbors(cur_node, colored)
            if colored and cur_node in orphaned_nodes:
                extra_neighbors = orphaned_nodes - set([cur_node])
                neighbors.extend(list(extra_neighbors))

            for other_node in neighbors:
                weight = self._edge_weight(cur_node, other_node, colored)

                if dist[other_node] > dist[cur_node] + weight:
                    dist[other_node] = dist[cur_node] + weight
                    parent[other_node] = cur_node
                    queue.insert((other_node, not colored), dist[other_node])

        if not found:
            return None

        path = [dst]
        cur_node = dst
        while cur_node != src:
            path.append(parent[cur_node])
            cur_node = parent[cur_node]
        return path[::-1]
