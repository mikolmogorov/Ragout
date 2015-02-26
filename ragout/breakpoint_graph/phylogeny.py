#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module solves "Half-breakpoint state parsimony"
problem
"""

import math
from collections import defaultdict
import logging

from ragout.parsers.phylogeny_parser import (parse_tree, PhyloException)
logger = logging.getLogger()

class Phylogeny:
    """
    Represents phylogenetic tree and scores it with
    given half-breakpoint states
    """
    def __init__(self, recipe):
        self.tree_string = recipe["tree"]
        self.tree = parse_tree(self.tree_string)
        self.scale_branches()

    def scale_branches(self):
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
        logger.debug("Branch lengths: {0}, mu = {1}".format(lengths, self.mu))

    def estimate_tree(self, leaf_states):
        """
        Scoring with DP (see algorithm description in the paper)
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


def _median(values):
    sorted_values = sorted(values)
    return sorted_values[(len(values) - 1) / 2]
