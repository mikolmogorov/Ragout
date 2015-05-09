#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module tries to detect chimerias in
input scaffolds
"""

from __future__ import print_function
import logging
from collections import defaultdict
from copy import copy, deepcopy

from ragout.breakpoint_graph.breakpoint_graph import BreakpointGraph

logger = logging.getLogger()

class ChimeraDetector(object):
    def __init__(self, breakpoint_graph, perm_container, target_seqs):
        self.graph = breakpoint_graph
        self.target_seqs = target_seqs
        self._get_chimeric_adj()
        self._get_contig_cuts(perm_container.target_perms)

    def _get_chimeric_adj(self):
        logger.info("Detecting chimeric sequences")
        chimeric_adj = set()

        subgraphs = self.graph.connected_components()
        for subgr in subgraphs:
            if len(subgr.bp_graph) > 100:
                logger.debug("Processing component of size {0}"
                             .format(len(subgr.bp_graph)))

            for (u, v, data) in subgr.bp_graph.edges_iter(data=True):
                genomes = subgr.supporting_genomes(u, v)
                if set(genomes) != set([self.graph.target]):
                    continue
                if subgr.alternating_cycle(u, v, False):
                    continue

                gap_seq = (self.target_seqs[data["chr_name"]]
                                           [data["start"]:data["end"]])
                ns_rate = (float(gap_seq.upper().count("N")) / len(gap_seq)
                           if len(gap_seq) else 0)
                if ns_rate < 0.1 and len(subgr.bp_graph) == 4:
                    if self._valid_2break(subgr.bp_graph, (u, v)):
                        logger.debug(ns_rate)
                        for node in subgr.bp_graph.nodes():
                            self.graph.add_debug_node(node)
                        continue
                chimeric_adj.add(tuple(sorted([u, v])))

        self.graph.debug_output()
        self.chimeric_adj = chimeric_adj

    def _valid_2break(self, bp_graph, red_edge):
        assert len(bp_graph) == 4
        red_1, red_2 = red_edge
        cand_1, cand_2 = tuple(set(bp_graph.nodes()) - set(red_edge))
        if abs(cand_1) == abs(cand_2):
            return False

        if bp_graph.has_edge(red_1, cand_1):
            known_1 = red_1, cand_1
            known_2 = red_2, cand_2
        else:
            known_1 = red_1, cand_2
            known_2 = red_2, cand_1

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

    def _get_contig_cuts(self, permutations):
        cuts = defaultdict(list)
        for perm in permutations:
            for pb, nb in perm.iter_pairs():
                adj = tuple(sorted([-pb.signed_id(), nb.signed_id()]))
                if adj in self.chimeric_adj:
                    cut = (nb.start + pb.end) / 2
                    cuts[perm.chr_name].append(cut)

        self.cuts = cuts

    def fix_container(self, perm_container):
        perm_container.target_perms = self._cut_permutations(self.cuts,
                                                perm_container.target_perms)
        perm_container.filter_indels(True)

    #TODO: refactoring
    def _cut_permutations(self, cuts, permutations):
        new_perms = []
        num_chim_perms = 0
        num_cuts = 0
        num_lost = 0
        for perm in permutations:
            if perm.chr_name not in cuts:
                new_perms.append(perm)
                continue

            #logger.debug("Original {0}".format(perm))
            perm_cuts = copy(cuts[perm.chr_name])
            perm_cuts.append(perm.chr_len)
            cur_perm = deepcopy(perm)
            cur_perm.blocks = []
            shift = 0

            num_chim_perms += 1
            num_cuts += len(perm_cuts) - 1

            for block in perm.blocks:
                if block.end <= perm_cuts[0]:
                    block.start -= shift
                    block.end -= shift
                    cur_perm.blocks.append(block)
                    continue

                if block.start < perm_cuts[0]:
                    num_lost += 1
                    continue

                #we have passed the current cut
                cur_perm.chr_name = (cur_perm.chr_name +
                                   "[{0}:{1}]".format(shift, perm_cuts[0]))
                cur_perm.chr_len = perm_cuts[0] - shift
                new_perms.append(cur_perm)
                #logger.debug(cur_perm)

                shift = perm_cuts[0]
                perm_cuts.pop(0)

                cur_perm = deepcopy(perm)
                block.start -= shift
                block.end -= shift
                cur_perm.blocks = [block]

            cur_perm.chr_name = (cur_perm.chr_name +
                               "[{0}:{1}]".format(shift, perm_cuts[0]))
            cur_perm.chr_len = perm_cuts[0] - shift
            new_perms.append(cur_perm)
            #logger.debug(cur_perm)

        logger.debug("Chimera Detector: {0} cuts made in {1} sequences"
                        .format(num_cuts, num_chim_perms))
        logger.debug("Lost {0} blocks".format(num_lost))
        return new_perms
