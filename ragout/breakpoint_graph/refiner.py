#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
Various algorithms for breakpoint graph processing
"""

from __future__ import print_function
import logging
from collections import defaultdict, namedtuple
from itertools import product

import networkx as nx

from ragout.shared.priority_queue import PriorityQueue

logger = logging.getLogger()
Adjacency = namedtuple("Adjacency", ["block", "distance", "supporting_genomes"])

class AdjacencyRefiner(object):
    def __init__(self, breakpoint_graph, phylogeny, perm_container):
        self.bp_graph = breakpoint_graph
        self.phylogeny = phylogeny
        self.perm_container = perm_container

    def refine_adjacencies(self, scaffolds):
        """
        Finding adjacencies consisten with previous iteration
        """
        orphaned_nodes = self.bp_graph.get_orphaned_nodes()
        weighted_graph = self.bp_graph.make_weighted(self.phylogeny)
        trusted_adj, mandatory_adj = self._get_trusted_adjacencies(scaffolds)
        chosen_edges = _get_path_cover(weighted_graph, trusted_adj,
                                       mandatory_adj, orphaned_nodes)

        adjacencies = {}
        for node_1, node_2 in chosen_edges:
            #infinity edges correspond to joined chromosome ends -- ignore them
            if self.bp_graph.is_infinity(node_1, node_2):
                continue

            distance = self.bp_graph.get_distance(node_1, node_2)
            supporting_genomes = self.bp_graph \
                                        .supporting_genomes(node_1, node_2)
            assert abs(node_1) != abs(node_2)
            adjacencies[node_1] = Adjacency(node_2, distance,
                                            supporting_genomes)
            adjacencies[node_2] = Adjacency(node_1, distance,
                                            supporting_genomes)
        return adjacencies

    def _get_trusted_adjacencies(self, prev_scaffolds):
        """
        Get trusted adjaencies from previous iteration
        """
        trusted_adj = []
        mandatory_adj = []
        perm_by_id = {perm.chr_name : perm for perm in
                      self.perm_container.target_perms}

        for scf in prev_scaffolds:
            for prev_cont, next_cont in zip(scf.contigs[:-1], scf.contigs[1:]):
                if (prev_cont.seq_name in perm_by_id and
                    next_cont.seq_name in perm_by_id):
                    left_blocks = perm_by_id[prev_cont.seq_name].blocks
                    left = (left_blocks[-1].signed_id() if prev_cont.sign > 0
                            else -left_blocks[0].signed_id())

                    right_blocks = perm_by_id[next_cont.seq_name].blocks
                    right = (right_blocks[0].signed_id() if next_cont.sign > 0
                             else -right_blocks[-1].signed_id())

                    #TODO: fix it
                    assert prev_cont.seq_name != next_cont.seq_name
                    trusted_adj.append((-left, right))

            for perm in self.perm_container.target_perms:
                mandatory_adj.append((-perm.blocks[0].signed_id(),
                                      perm.blocks[-1].signed_id()))

        return trusted_adj, mandatory_adj


def _get_path_cover(graph, trusted_adj, mandatory_adj, candidate_nodes):
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
