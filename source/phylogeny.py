import math
from Bio import Phylo
from cStringIO import StringIO
from collections import defaultdict


class Phylogeny:
    def __init__(self, newick_string):
        self.tree = Phylo.read(StringIO(newick_string), "newick")
        self.validate_tree()

    #TODO
    def validate_tree(self):
        pass

    def estimate_tree(self, adjacencies):
        dp = tree_score(self.tree, adjacencies)
        return dp, None, None


def tree_score(tree, leaf_states):
    all_states = set(leaf_states.values())

    def branch_score(root, child, branch):
        MU = 1
        if root == child:
            return 0.0
        else:
            length = max(branch, 0.0000001)
            return math.exp(-MU * length)

    def rec_helper(root):
        if root.is_terminal():
            leaf_score = lambda s: 0.0 if s == leaf_states[root.name] else float("inf")
            return {s : leaf_score(s) for s in all_states}

        nodes_scores = {}
        for node in root.clades:
            nodes_scores[node] = rec_helper(node)

        root_scores = defaultdict(float)
        for root_state in all_states:
            for node in root.clades:
                min_score = float("inf")
                for child_state in all_states:
                    score = (nodes_scores[node][child_state] +
                            branch_score(root_state, child_state, node.branch_length))
                    min_score = min(min_score, score)
                root_scores[root_state] += min_score

        return root_scores

    return min(rec_helper(tree.clade).values())
