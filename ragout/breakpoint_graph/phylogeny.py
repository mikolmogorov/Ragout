#This module solves "Half-breakpoint state parsimony"
#problem
#####################################################

import math
from collections import defaultdict
try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO

from Bio import Phylo

#PUBLIC:
####################################################

class PhyloException(Exception):
    pass

#Represents phylogenetic tree and scores it with 
#given half-breakpoint states
class Phylogeny:
    def __init__(self, recipe):
        self.tree = Phylo.read(StringIO(recipe.tree), "newick")
        self.tree.clade.branch_length = 0
        self.tree_string = recipe.tree
        genomes = dict(recipe.references.items() + recipe.targets.items())
        self.validate_tree(genomes)

    def validate_tree(self, recipe_genomes):
        tree_genomes = set(map(lambda n: n.name, self.tree.get_terminals()))
        if tree_genomes != set(recipe_genomes.keys()):
            raise PhyloException("Some genomes are missing "
                                 "from the tree/description")

    def estimate_tree(self, adjacencies):
        return _tree_score(self.tree, adjacencies)

    #PRIVATE:
####################################################

#scoring with DP (see algorithm description in the paper)
def _tree_score(tree, leaf_states):
    all_states = set(leaf_states.values())

    #score of a tree branch
    def branch_score(root, child, branch):
        MU = 1
        if root == child:
            return 0.0
        else:
            length = max(branch, 0.0000001)
            return math.exp(-MU * length)

    #recursive
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
