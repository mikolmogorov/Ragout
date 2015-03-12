#(c) 2013-2015 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module infers phylogenetic tree based on
breakpoints data
"""

from collections import defaultdict
from itertools import combinations, product

from newick.tree import Leaf, Tree

class TreeBuilder:
    def __init__(self, permutations):
        self.perms_by_genome = defaultdict(list)
        for perm in permutations:
            self.perms_by_genome[perm.genome_name].append(perm)

    def _get_distance(self, genome_1, genome_2):
        breakpoints_1 = set()
        n_blocks_1 = 0
        for perm in self.perms_by_genome[genome_1]:
            n_blocks_1 += len(perm.blocks)
            for bl_1, bl_2 in zip(perm.blocks[:-1], perm.blocks[1:]):
                bp = sorted([-bl_1.signed_id(), bl_2.signed_id()])
                breakpoints_1.add(tuple(bp))

        breakpoints_2 = set()
        n_blocks_2 = 0
        for perm in self.perms_by_genome[genome_2]:
            n_blocks_2 += len(perm.blocks)
            for bl_1, bl_2 in zip(perm.blocks[:-1], perm.blocks[1:]):
                bp = sorted([-bl_1.signed_id(), bl_2.signed_id()])
                breakpoints_2.add(tuple(bp))

        #print("breakpoints", len(breakpoints_1), len(breakpoints_2))
        #print("Differences", breakpoints_1 ^ breakpoints_2)
        #return (max(n_blocks_1, n_blocks_2) -
        #        len(breakpoints_1 & breakpoints_2) - 1)
        return (min(len(breakpoints_1), len(breakpoints_2)) -
                len(breakpoints_1 & breakpoints_2))


    def build(self):
        genomes = self.perms_by_genome.keys()
        distances = defaultdict(lambda : {})
        for g_1, g_2 in combinations(genomes, 2):
            distances[g_1][g_2] = self._get_distance(g_1, g_2)
            distances[g_2][g_1] = distances[g_1][g_2]

        def tree_dist(tree_1, tree_2):
            leaves_1 = tree_1.get_leaves_identifiers()
            leaves_2 = tree_2.get_leaves_identifiers()
            total = 0
            for l_1, l_2 in product(leaves_1, leaves_2):
                total += distances[l_1][l_2]
            return float(total) / len(leaves_1) / len(leaves_2)

        trees_left = set(map(Leaf, genomes))
        while len(trees_left) > 1:
            #determine two closest ones
            lowest_dst = float("inf")
            lowest_pair = None
            for t_1, t_2 in combinations(trees_left, 2):
                dst = tree_dist(t_1, t_2)
                if dst < lowest_dst:
                    lowest_dst = dst
                    lowest_pair = (t_1, t_2)

            new_tree = Tree()
            new_tree.add_edge((lowest_pair[0], None, lowest_dst / 2))
            new_tree.add_edge((lowest_pair[1], None, lowest_dst / 2))
            trees_left.add(new_tree)
            trees_left.remove(lowest_pair[0])
            trees_left.remove(lowest_pair[1])

        tree = list(trees_left)[0]
        print(tree)
