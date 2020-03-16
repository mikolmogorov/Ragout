#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module tries to detect missassembled adjacencies
in input sequences and breaks them if neccesary
"""

from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
import logging
from collections import defaultdict, namedtuple
from copy import copy, deepcopy
from ragout.six.moves import range


logger = logging.getLogger()
ContigBreak = namedtuple("ContigBreak", ["seq_name", "begin", "end", "good"])

class ChimeraDetector(object):
    def __init__(self, breakpoint_graphs, run_stages, target_seqs):
        logger.info("Detecting chimeric adjacencies")
        self.bp_graphs = breakpoint_graphs
        self.run_stages = run_stages
        self.target_seqs = target_seqs
        self._make_hierarchical_breaks()

    def _make_hierarchical_breaks(self):
        """
        Determines where and at what synteny blocks scale to break each contig
        """
        seq_cuts = defaultdict(lambda : defaultdict(list))

        #extracting and grouping by sequence
        for stage in self.run_stages:
            logger.debug(">>With block size: %d", stage.block_size)
            breaks = self._get_contig_breaks(self.bp_graphs[stage])
            for br in breaks:
                seq_cuts[br.seq_name][stage].append(br)

        hierarchical_cuts = defaultdict(lambda : defaultdict(list))
        for seq_name in seq_cuts:
            for i in range(len(self.run_stages)):
                top_stage = self.run_stages[i]

                for top_break in seq_cuts[seq_name][top_stage]:
                    if top_break.good:
                        continue

                    adjusted_break = (top_break.begin, top_break.end)
                    for j in range(i + 1, len(self.run_stages)):
                        lower_stage = self.run_stages[j]
                        #check if there is overlapping cut
                        for lower_break in seq_cuts[seq_name][lower_stage]:
                            ovlp_left = max(adjusted_break[0], lower_break.begin)
                            ovlp_right = min(adjusted_break[1], lower_break.end)
                            #if so, update current and go down
                            if ovlp_right >= ovlp_left:
                                adjusted_break = (ovlp_left, ovlp_right)
                                break

                    break_pos = self._optimal_break(seq_name, *adjusted_break)
                    hierarchical_cuts[seq_name][top_stage].append(break_pos)
        self.hierarchical_cuts = hierarchical_cuts

    def _optimal_break(self, seq_name, break_start, break_end):
        """
        Finds the longet run of Ns within given range
        """
        seq = self.target_seqs[seq_name]
        cur_pos = break_start
        cur_len = 0
        max_pos = break_start
        max_len = 0
        for i in range(break_start, break_end):
            if seq[i].upper() == "N":
                cur_len += 1
            if seq[i].upper() != "N" or i == break_end - 1:
                if max_len < cur_len:
                    max_len = cur_len
                    max_pos = cur_pos
                cur_pos = i + 1
                cur_len = 0
        return max_pos + max_len // 2

    def _get_contig_breaks(self, bp_graph):
        """
        Detects chimeric adjacencies
        """
        seq_cuts = []

        subgraphs = bp_graph.connected_components()
        num_target_adj = 0
        num_unique_adj = 0
        num_removed_adj = 0
        for subgr in subgraphs:
            if len(subgr.bp_graph) > 100:
                logger.debug("Processing component of size %d",
                             len(subgr.bp_graph))

            for (u, v, data) in subgr.bp_graph.edges(data=True):
                if data["genome_id"] != subgr.target:
                    continue

                num_target_adj += 1
                genomes = set(subgr.genomes_support(u, v))
                if (genomes != set([bp_graph.target])
                    and not subgr.is_infinity(u, v)):
                    continue

                num_unique_adj += 1
                seq_name, start, end = data["chr_name"], data["start"], data["end"]
                if subgr.alternating_cycle(u, v) is not None:
                    seq_cuts.append(ContigBreak(seq_name, start, end, True))
                    continue

                seq_cuts.append(ContigBreak(seq_name, start, end, False))
                num_removed_adj += 1

        #bp_graph.debug_output()
        logger.debug("Checking %d target adjacencies", num_target_adj)
        logger.debug("Found %d target-specific adjacencies, "
                     "%d broken as chimeric", num_unique_adj, num_removed_adj)
        return seq_cuts

    def _valid_2break(self, bp_graph, red_edge):
        """
        Checks if there is a valid 2-break through the given red edge
        """
        assert len(bp_graph) == 4
        red_1, red_2 = red_edge
        cand_1, cand_2 = tuple(set(bp_graph.nodes) - set(red_edge))
        if abs(cand_1) == abs(cand_2):
            return False

        if bp_graph.has_edge(red_1, cand_1):
            if not bp_graph.has_edge(red_2, cand_2):
                return False
            known_1 = red_1, cand_1
            known_2 = red_2, cand_2
        elif bp_graph.has_edge(red_1, cand_2):
            if not bp_graph.has_edge(red_2, cand_1):
                return False
            known_1 = red_1, cand_2
            known_2 = red_2, cand_1
        else:
            return False

        chr_1 = {}
        for data in bp_graph[known_1[0]][known_1[1]].values():
            chr_1[data["genome_id"]] = data["chr_name"]
        chr_2 = {}
        for data in bp_graph[known_2[0]][known_2[1]].values():
            chr_2[data["genome_id"]] = data["chr_name"]
        common_genomes = set(chr_1.keys()).intersection(list(chr_2.keys()))
        for genome in common_genomes:
            if chr_1[genome] != chr_2[genome]:
                return False

        return True

    def break_contigs(self, perm_container, block_sizes):
        """
        Breaks contigs in inferred cut positions
        """
        logger.info("Removing chimeric adjacencies")

        new_container = deepcopy(perm_container)
        new_target_perms = []
        num_breaks = 0
        num_chimeras = 0

        for perm in new_container.target_perms:
            break_points = []
            for size in block_sizes:
                break_points.extend(self.hierarchical_cuts[perm.chr_name][size])
            break_points = list(set(break_points))
            num_breaks += len(break_points)
            if not break_points:
                new_target_perms.append(perm)
            else:
                num_chimeras += 1
                new_target_perms.extend(_break_permutation(perm, break_points))

        logger.debug("Chimera Detector: %d cuts made in %d sequences",
                     num_breaks, num_chimeras)
        new_container.target_perms = new_target_perms
        return new_container


def _break_permutation(permutation, break_points):
    broken_perms = []

    cuts_stack = copy(sorted(break_points))
    cuts_stack.append(permutation.seq_len)
    current_perm = deepcopy(permutation)
    current_perm.blocks = []
    shift = 0

    for block in permutation.blocks:
        if block.end <= cuts_stack[0]:
            block.start -= shift
            block.end -= shift
            current_perm.blocks.append(block)
            continue

        if block.start < cuts_stack[0]:
            block.start = cuts_stack[0]

        #we have passed the current cut
        current_perm.seq_start = shift
        current_perm.seq_end = cuts_stack[0]
        if current_perm.blocks:
            broken_perms.append(current_perm)

        shift = cuts_stack[0]
        cuts_stack.pop(0)

        current_perm = deepcopy(permutation)
        block.start -= shift
        block.end -= shift
        current_perm.blocks = [block]

    current_perm.seq_start = shift
    current_perm.seq_end = cuts_stack[0]
    if current_perm.blocks:
        broken_perms.append(current_perm)

    return broken_perms
