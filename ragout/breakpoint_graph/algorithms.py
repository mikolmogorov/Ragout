#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
Various algorithms for breakpoint graph processing
"""

from __future__ import print_function
import heapq
import logging
from collections import defaultdict
from itertools import product

import networkx as nx

logger = logging.getLogger()


def min_weight_matching(graph):
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


def alternating_cycle(graph, node_1, node_2, target_id):
    """
    Determines if there is a cycle of alternating colors
    that goes through the given edge
    """
    def get_genome_ids((u, v)):
        return list(map(lambda e: e["genome_id"], graph[u][v].values()))

    simple_graph = nx.Graph()
    for (u, v) in graph.edges_iter():
        simple_graph.add_edge(u, v)
    assert not simple_graph.has_edge(node_1, node_2)

    good_path = False
    for path in nx.all_simple_paths(simple_graph, node_1, node_2):
        if len(path) % 2 == 1:
            continue

        edges = list(zip(path[:-1], path[1:]))
        odd_colors = list(map(get_genome_ids, edges[0::2]))
        even_colors = list(map(get_genome_ids, edges[1::2]))

        if not all(map(lambda e: set(e) == set([target_id]), even_colors)):
            continue

        common_genomes = set(odd_colors[0])
        for edge_colors in odd_colors:
            common_genomes = common_genomes.intersection(edge_colors)

        if common_genomes:
            good_path = True
            break

    return len(path) / 2 if good_path else None


def get_orphaned_nodes(graph, target_id):
    """
    Get nodes suspected to rearrangements
    """
    logger.debug("Getting candidate nodes")
    #graph = graph.copy()
    candidate_nodes = set()

    subgraphs = list(nx.connected_component_subgraphs(graph))
    for subgr in subgraphs:
        known_nodes = set(subgr.nodes())
        for v1, v2, data in subgr.edges_iter(data=True):
            genome_ids = list(map(lambda e: e["genome_id"],
                                  subgr[v1][v2].values()))
            if target_id in genome_ids:
                known_nodes.discard(v1)
                known_nodes.discard(v2)

        if len(known_nodes) == 2:
            node_1, node_2 = tuple(known_nodes)
            if subgr.has_edge(node_1, node_2):
                continue

            cycle = alternating_cycle(subgr, node_1, node_2, target_id)
            if (abs(node_1) != abs(node_2) and cycle is not None):
                candidate_nodes.add(node_1)
                candidate_nodes.add(node_2)

    return candidate_nodes

    """
    logger.debug("Modifying graph")
    counter = 0
    for n1, n2 in product(candidate_nodes, repeat=2):
        if abs(n1) != abs(n2):
            graph.add_edge(n1, n2, candidate=True)
            counter += 1
    logger.debug("Added {0} extra edges".format(counter))

    return graph
    """


def get_path_cover(graph, trusted_adj, mandatory_adj, candidate_nodes):
    logger.debug("Computing path cover")
    adjacencies = []
    prohibited_nodes = set()
    black_adj = {}
    for (u, v) in mandatory_adj:
        black_adj[u] = v
        black_adj[v] = u

    for adj in trusted_adj:
        prohibited_nodes.add(abs(adj[0]))
        prohibited_nodes.add(abs(adj[1]))

    for (adj_left, adj_right) in trusted_adj:
        p = _shortest_path(graph, adj_left, adj_right, prohibited_nodes,
                           black_adj, candidate_nodes)
        #logger.debug(p)
        if not p:
            p = [adj_left, adj_right]

        #assert len(p) % 2 == 0
        if len(p) % 2 != 0:
            print(p)
            p = [adj_left, adj_right]

        for i in xrange(len(p) / 2):
            adj_left, adj_right = p[i * 2], p[i * 2 + 1]
            adjacencies.append((adj_left, adj_right))
            prohibited_nodes.add(abs(adj_left))
            prohibited_nodes.add(abs(adj_right))

    return adjacencies


def _shortest_path(graph, src, dst, prohibited_nodes,
                   black_adj, candidate_nodes):
    """
    Finds shortest path wrt to restricted nodes
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

        if cur_node != src and abs(cur_node) in prohibited_nodes:
            continue
        if not colored and -cur_node not in black_adj:
            continue

        neighbors = (graph.neighbors(cur_node) if colored
                     else [-black_adj[-cur_node]])
        if colored and cur_node in candidate_nodes:
            extra_neighbors = candidate_nodes - set([cur_node])
            neighbors.extend(list(extra_neighbors))

        for other_node in neighbors:
            if colored:
                if graph.has_edge(cur_node, other_node):
                    weight = graph[cur_node][other_node]["weight"]
                else:
                    weight = 1
            else:
                weight = 0

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


class PriorityQueue(object):
    """
    Priority queue based on heap, capable of inserting a new node with
    desired priority, updating the priority of an existing node and deleting
    an abitrary node while keeping invariant
    """

    def __init__(self):
        """
        if 'heap' is not empty, make sure it's heapified
        """
        self.heap = []
        self.entry_finder = {}
        self.REMOVED = "<remove_marker>"
        self.length = 0

    def get_length(self):
        return self.length

    def insert(self, node, priority=0):
        """
        'entry_finder' bookkeeps all valid entries, which are bonded in
        'heap'. Changing an entry in either leads to changes in both.
        """
        if node in self.entry_finder:
            self.delete(node)
        entry = [priority, node]
        self.entry_finder[node] = entry
        heapq.heappush(self.heap, entry)
        self.length += 1

    def delete(self, node):
        """
        Instead of breaking invariant by direct removal of an entry, mark
        the entry as "REMOVED" in 'heap' and remove it from 'entry_finder'.
        Logic in 'pop()' properly takes care of the deleted nodes.
        """
        entry = self.entry_finder.pop(node)
        entry[-1] = self.REMOVED
        self.length -= 1
        return entry[0]

    def pop(self):
        """
        Any popped node marked by "REMOVED" does not return, the deleted
        nodes might be popped or still in heap, either case is fine.
        """
        while self.heap:
            priority, node = heapq.heappop(self.heap)
            if node is not self.REMOVED:
                del self.entry_finder[node]
                self.length -= 1
                return priority, node
        raise KeyError('pop from an empty priority queue')
