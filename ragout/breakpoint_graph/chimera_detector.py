#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module tries to detect missassembled adjacencies
in input sequences and breaks them if neccesary
"""

from __future__ import print_function
import logging
from collections import defaultdict, namedtuple
from copy import copy, deepcopy

from ragout.breakpoint_graph.breakpoint_graph import BreakpointGraph

logger = logging.getLogger()
ContigBreak = namedtuple("ContigBreak", ["seq_name", "begin", "end", "good"])

class ChimeraDetector(object):
    def __init__(self, breakpoint_graphs, run_stages):
        logger.debug("Detecting chimeric adjacencies")
        self.bp_graphs = breakpoint_graphs
        self.run_stages = run_stages
        self._make_hierarchical_breaks()

    def _make_hierarchical_breaks(self):
        """
        Determines where and at what synteny blocks scale to break each contig
        """
        seq_cuts = defaultdict(lambda : defaultdict(list))

        #extracting and grouping by sequence
        for stage in self.run_stages:
            breaks = self._get_contig_breaks(self.bp_graphs[stage])
            for br in breaks:
                seq_cuts[br.seq_name][stage].append(br)

        hierarchical_cuts = defaultdict(lambda : defaultdict(list))
        for seq_name in seq_cuts:
            for i in xrange(len(self.run_stages)):
                top_stage = self.run_stages[i]

                for top_break in seq_cuts[seq_name][top_stage]:
                    if top_break.good:
                        continue

                    adjusted_break = (top_break.begin, top_break.end)
                    for j in xrange(i + 1, len(self.run_stages)):
                        lower_stage = self.run_stages[j]
                        #check if there is overlapping cut
                        for lower_break in seq_cuts[seq_name][lower_stage]:
                            ovlp_left = max(adjusted_break[0], lower_break.begin)
                            ovlp_right = min(adjusted_break[1], lower_break.end)
                            #if so, update current and go down
                            if ovlp_right >= ovlp_left:
                                adjusted_break = (ovlp_left, ovlp_right)
                                break

                    break_pos = adjusted_break[0]
                    hierarchical_cuts[seq_name][top_stage].append(break_pos)
        self.hierarchical_cuts = hierarchical_cuts

    def _get_contig_breaks(self, bp_graph):
        """
        Detects chimeric adjacencies
        """
        seq_cuts = []

        subgraphs = bp_graph.connected_components()
        for subgr in subgraphs:
            if len(subgr.bp_graph) > 100:
                logger.debug("Processing component of size {0}"
                             .format(len(subgr.bp_graph)))

            for (u, v, data) in subgr.bp_graph.edges_iter(data=True):
                if data["genome_id"] != subgr.target:
                    continue

                genomes = set(subgr.genomes_support(u, v))
                if (genomes != set([bp_graph.target])
                    and not subgr.is_infinity(u, v)):
                    continue

                seq_name, start, end = data["chr_name"], data["start"], data["end"]
                if subgr.alternating_cycle(u, v) is not None:
                    seq_cuts.append(ContigBreak(seq_name, start, end, True))
                    continue

                seq_cuts.append(ContigBreak(seq_name, start, end, False))

        #bp_graph.debug_output()
        return seq_cuts

    def _valid_2break(self, bp_graph, red_edge):
        """
        Checks if there is a valid 2-break through the given red edge
        """
        assert len(bp_graph) == 4
        red_1, red_2 = red_edge
        cand_1, cand_2 = tuple(set(bp_graph.nodes()) - set(red_edge))
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
        common_genomes = set(chr_1.keys()).intersection(chr_2.keys())
        for genome in common_genomes:
            if chr_1[genome] != chr_2[genome]:
                return False

        return True

    def break_contigs(self, perm_container, block_sizes):
        """
        Breaks contigs in inferred cut positions
        """
        new_container = deepcopy(perm_container)
        new_container.target_perms = \
                self._cut_permutations(new_container.target_perms, block_sizes)
        return new_container

    def _cut_permutations(self, permutations, block_sizes):
        """
        Actually breaks these contigs
        """
        new_perms = []
        num_chim_perms = 0
        num_cuts = 0
        num_lost = 0
        for perm in permutations:
            cuts = []
            for size in block_sizes:
                cuts.extend(self.hierarchical_cuts[perm.chr_name][size])
            cuts = list(set(cuts))
            if not cuts:
                new_perms.append(perm)
                continue

            #logger.debug("Original {0}".format(perm))
            cuts_stack = copy(sorted(cuts))
            cuts_stack.append(perm.seq_len)
            cur_perm = deepcopy(perm)
            cur_perm.blocks = []
            shift = 0

            num_chim_perms += 1
            num_cuts += len(cuts_stack) - 1

            for block in perm.blocks:
                if block.end <= cuts_stack[0]:
                    block.start -= shift
                    block.end -= shift
                    cur_perm.blocks.append(block)
                    continue

                if block.start < cuts_stack[0]:
                    num_lost += 1
                    continue

                #we have passed the current cut
                cur_perm.seq_start = shift
                cur_perm.seq_end = cuts_stack[0]
                if cur_perm.blocks:
                    new_perms.append(cur_perm)
                #logger.debug(cur_perm)

                shift = cuts_stack[0]
                cuts_stack.pop(0)

                cur_perm = deepcopy(perm)
                block.start -= shift
                block.end -= shift
                cur_perm.blocks = [block]

            cur_perm.seq_start = shift
            cur_perm.seq_end = cuts_stack[0]
            if cur_perm.blocks:
                new_perms.append(cur_perm)
            #logger.debug(cur_perm)

        logger.debug("Chimera Detector: {0} cuts made in {1} sequences"
                        .format(num_cuts, num_chim_perms))
        logger.debug("Lost {0} blocks".format(num_lost))
        return new_perms
