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
        #orphaned_nodes = self.bp_graph.get_orphaned_nodes()
        adj_graph = AdjacencyGraph(self.bp_graph, self.phylogeny)
        trusted_adj = self._get_trusted_adjacencies(scaffolds)
        return self._path_cover(adj_graph, trusted_adj, set())

    def _get_trusted_adjacencies(self, prev_scaffolds):
        """
        Get trusted adjaencies from previous iteration
        """
        trusted_adj = []
        for scf in prev_scaffolds:
            for prev_cont, next_cont in zip(scf.contigs[:-1], scf.contigs[1:]):
                left, right = prev_cont.right_end(), next_cont.left_end()
                flank = prev_cont.right_gap() + next_cont.left_gap()
                distance = prev_cont.link.gap + flank
                genomes = prev_cont.link.supporting_genomes

                trusted_adj.append((left, right, distance, genomes))

        return trusted_adj

    def _path_cover(self, adj_graph, trusted_adj, orphaned_nodes):
        logger.debug("Computing path cover")
        adjacencies = {}
        prohibited_nodes = set()
        for adj in trusted_adj:
            prohibited_nodes.add(adj[0])
            prohibited_nodes.add(adj[1])

        def add_adj(left, right, dist, genomes):
            prohibited_nodes.add(left)
            prohibited_nodes.add(right)

            adjacencies[left] = Adjacency(adj_right, dist, genomes)
            adjacencies[right] = Adjacency(adj_left, dist, genomes)

        generated_adj = {}
        for (adj_left, adj_right, dist, genomes) in trusted_adj:
            p = adj_graph.shortest_path(adj_left, adj_right, prohibited_nodes,
                                        orphaned_nodes)
            logger.debug(p)
            if not p or len(p) % 2 == 1:
                #assert len(p) % 2 == 0
                add_adj(adj_left, adj_right, dist, genomes)
            else:
                for i in xrange(len(p) / 2):
                    adj_left, adj_right = p[i * 2], p[i * 2 + 1]
                    #assert abs(adj_left) != abs(adj_right)
                    if self.bp_graph.is_infinity(adj_left, adj_right):
                        continue
                    if abs(adj_left) == abs(adj_right):
                        continue

                    dist = self.bp_graph.get_distance(adj_left, adj_right)
                    genomes = self.bp_graph.supporting_genomes(adj_left,
                                                               adj_right)
                    add_adj(adj_left, adj_right, dist, genomes)

        return adjacencies
