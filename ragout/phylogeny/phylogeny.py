#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module provides some phylogeny-related functions
(includeing one for solving half-breakpoint state parsimony problem)
"""

from __future__ import absolute_import
from __future__ import division
import math
from collections import defaultdict
from itertools import chain
import logging

import networkx as nx

from ragout.parsers.phylogeny_parser import (parse_tree)
from ragout.phylogeny.inferer import TreeInferer
logger = logging.getLogger()

class Phylogeny:
    """
    Represents phylogenetic tree and scores it with
    given half-breakpoint states
    """
    def __init__(self, tree):
        self.tree = tree
        self.tree_string = str(tree)
        self._scale_branches()

    @classmethod
    def from_newick(phylo, newick_str):
        return phylo(parse_tree(newick_str))

    @classmethod
    def from_permutations(phylo, perm_container):
        ti = TreeInferer(perm_container)
        return phylo(ti.build())

    def _scale_branches(self):
        """
        Fits mu coefficient according to branch lengths
        to avoid underflows/overflows
        """
        def tree_length(node):
            if node.terminal:
                return []

            lengths = []
            for node, _bootstrap, length in node.get_edges():
                assert length is not None
                lengths.append(length)
                lengths.extend(tree_length(node))

            return lengths

        lengths = tree_length(self.tree)
        assert len(lengths)
        self.mu = float(1) / _median(lengths)
        logger.debug("Branch lengths: %s, mu = %f", str(lengths), self.mu)

    def estimate_tree(self, leaf_states):
        """
        Scores the tree with weighted parsimony procedure
        """
        all_states = set(leaf_states.values())

        #score of a tree branch
        def branch_score(parent, child, branch):
            if parent == child or child is None:
                return 0.0
            else:
                #prevent underflow
                length = max(branch, 0.0000001)
                #adding one to counter possibly small exp value 
                return 1.0 + math.exp(-self.mu * length)

        #recursive
        def rec_helper(root):
            if root.terminal:
                leaf_score = (lambda s: 0.0 if s == leaf_states[root.identifier]
                                            else float("inf"))
                return {s : leaf_score(s) for s in all_states}

            nodes_scores = {}
            for node, _bootstrap, _length  in root.get_edges():
                nodes_scores[node] = rec_helper(node)

            root_scores = defaultdict(float)
            for root_state in all_states:
                for node, _bootstrap, branch_length in root.edges:
                    min_score = float("inf")
                    for child_state in all_states:
                        score = (nodes_scores[node][child_state] +
                                branch_score(root_state, child_state,
                                             branch_length))
                        min_score = min(min_score, score)
                    root_scores[root_state] += min_score

            return root_scores

        return min(rec_helper(self.tree).values())

    def terminals_dfs_order(self):
        """
        Returns terminal nodes' names in dfs order
        """
        def get_labels(root):
            if root.terminal:
                return [root.identifier]

            edges = sorted(root.get_edges(), key=lambda e: e[2], reverse=True)
            edges = sorted(edges, key=lambda e: e[0].terminal)
            #return list(chain(*map(lambda e: get_labels(e[0]), edges)))
            return list(chain(*[get_labels(e[0]) for e in edges]))

        return get_labels(self.tree)

    def leaves_by_distance(self, genome):
        """
        Returns leaves names sorted by the distance from
        the given genome.
        """
        graph = nx.Graph()
        start = [None]
        def rec_helper(root):
            if root.terminal and root.identifier == genome:
                start[0] = root
            if root.terminal:
                return
            for node, _bootstrap, branch_length in root.edges:
                graph.add_edge(root, node, weight=branch_length)
                rec_helper(node)

        rec_helper(self.tree)
        distances = nx.single_source_dijkstra_path_length(graph, start[0])
        leaves = [g for g in distances.keys()
                  if g.terminal and g.identifier != genome]

        #return list(map(lambda n: n.identifier,
        #                sorted(leaves, key=distances.get)))
        return [n.identifier for n in
                sorted(leaves, key=distances.get)]


def _median(values):
    sorted_values = sorted(values)
    return sorted_values[(len(values) - 1) // 2]
