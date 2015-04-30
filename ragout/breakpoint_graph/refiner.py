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

from ragout.breakpoint_graph.adjacency_graph import AdjacencyGraph

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
        adj_graph = AdjacencyGraph(self.bp_graph, self.phylogeny)
        trusted_adj = self._get_trusted_adjacencies(scaffolds)
        chosen_edges = _get_path_cover(adj_graph, trusted_adj, orphaned_nodes)

        adjacencies = {}
        for node_1, node_2 in chosen_edges:
            #infinity edges correspond to joined chromosome ends -- ignore them
            if self.bp_graph.is_infinity(node_1, node_2):
                continue
            #TODO: fix this
            if abs(node_1) == abs(node_2):
                continue

            distance = self.bp_graph.get_distance(node_1, node_2)
            supporting_genomes = self.bp_graph \
                                        .supporting_genomes(node_1, node_2)
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

                    assert prev_cont.seq_name != next_cont.seq_name
                    trusted_adj.append((-left, right))

        return trusted_adj


def _get_path_cover(adj_graph, trusted_adj, orphaned_nodes):
    logger.debug("Computing path cover")
    adjacencies = []
    prohibited_nodes = set()
    for adj in trusted_adj:
        prohibited_nodes.add(adj[0])
        prohibited_nodes.add(adj[1])

    for (adj_left, adj_right) in trusted_adj:
        p = adj_graph.shortest_path(adj_left, adj_right, prohibited_nodes,
                                    orphaned_nodes)
        #logger.debug(p)
        #TODO: fix %2
        if not p or len(p) % 2 == 1:
            p = [adj_left, adj_right]

        for i in xrange(len(p) / 2):
            adj_left, adj_right = p[i * 2], p[i * 2 + 1]
            adjacencies.append((adj_left, adj_right))
            prohibited_nodes.add(adj_left)
            prohibited_nodes.add(adj_right)

    return adjacencies
