#(c) 2013-2015 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module infers phylogenetic tree based on
breakpoints data
"""

from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from collections import defaultdict
from itertools import (combinations,
                       combinations_with_replacement, chain)

from ragout.newick.tree import Leaf, Tree
from ragout.six.moves import map
from ragout.six.moves import zip

class TreeInferer:
    def __init__(self, perm_container):
        self.perms_by_genome = defaultdict(list)
        for perm in chain(perm_container.ref_perms,
                          perm_container.target_perms):
            self.perms_by_genome[perm.genome_name].append(perm)

    def _genome_distance(self, genome_1, genome_2):
        """
        Calculates breakpoint distance between two genomes
        """
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

        return (min(len(breakpoints_1), len(breakpoints_2)) -
                len(breakpoints_1 & breakpoints_2))
        #return (max(n_blocks_1, n_blocks_2) - 
        #        len(breakpoints_1 & breakpoints_2) - 2)

    def build(self):
        """
        Implementation of neighbor-joining algorithm
        """
        MIN_LEN = 0.000001
        genomes = list(self.perms_by_genome.keys())
        taxas = list(map(Leaf, sorted(genomes)))
        for t in taxas:
            t.terminal = True

        distances = defaultdict(lambda : {})
        for t_1, t_2 in combinations_with_replacement(taxas, 2):
            distances[t_1][t_2] = self._genome_distance(t_1.identifier,
                                                        t_2.identifier)
            distances[t_2][t_1] = distances[t_1][t_2]

        def calc_q(taxas):
            q_matrix = defaultdict(lambda : {})
            for t_1, t_2 in combinations(taxas, 2):
                other_dist = 0
                for other_t in taxas:
                    other_dist += distances[t_1][other_t]
                    other_dist += distances[t_2][other_t]
                q_matrix[t_1][t_2] = ((len(taxas) - 2) * distances[t_1][t_2] -
                                     other_dist)
                q_matrix[t_2][t_1] = q_matrix[t_1][t_2]
            return q_matrix

        while len(taxas) > 1:
            #determine two closest ones
            q_matrix = calc_q(taxas)
            lowest_dst = float("inf")
            lowest_pair = None
            for t_1, t_2 in sorted(combinations(taxas, 2)):
                if q_matrix[t_1][t_2] < lowest_dst:
                    lowest_dst = q_matrix[t_1][t_2]
                    lowest_pair = (t_1, t_2)

            #calculate distances to new internal node from joined taxas
            new_taxa = Tree()
            new_taxa.terminal = False

            old_1, old_2 = sorted(lowest_pair)
            other_dist = 0
            for other_taxa in taxas:
                other_dist += distances[old_1][other_taxa]
                other_dist -= distances[old_2][other_taxa]
            div_dist = (0.5 / (len(taxas) - 2) * other_dist
                        if len(taxas) > 2 else 0)
            dist_1 = 0.5 * distances[old_1][old_2] + div_dist
            dist_2 = distances[old_1][old_2] - dist_1
            dist_1, dist_2 = max(MIN_LEN, dist_1), max(MIN_LEN, dist_2)

            new_taxa.add_edge((old_1, None, dist_1))
            new_taxa.add_edge((old_2, None, dist_2))
            taxas.remove(old_1)
            taxas.remove(old_2)

            for other_taxa in taxas:
                distances[new_taxa][other_taxa] = \
                    0.5 * (distances[old_1][other_taxa] +
                           distances[old_2][other_taxa] -
                           distances[old_1][old_2])
                distances[other_taxa][new_taxa] = distances[new_taxa][other_taxa]
            distances[new_taxa][new_taxa] = 0
            taxas.append(new_taxa)

        tree = list(taxas)[0]
        return tree
