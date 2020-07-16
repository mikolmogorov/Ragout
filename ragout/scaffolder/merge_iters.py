#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module provides some functions for
moving between two consecutive iterations
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from collections import defaultdict
from itertools import chain
import os
import logging
from copy import deepcopy, copy

import networkx as nx

from ragout.shared.debug import DebugConfig
from ragout.shared.datatypes import (Contig, Scaffold, Link,
                                     output_scaffolds_premutations)
from ragout.scaffolder.output_generator import output_links
from ragout.scaffolder.scaffolder import build_scaffolds
from ragout.breakpoint_graph.inferer import Adjacency
from ragout.breakpoint_graph.breakpoint_graph import GenChrPair
from ragout.six.moves import range
from ragout.six.moves import zip


logger = logging.getLogger()
debugger = DebugConfig.get_instance()


def merge_scaffolds(big_scaffolds, small_scaffolds, perm_container, rearrange):
    """
    Merges scaffold sets from different iterations. If rearrangements are allowed,
    tries to keep some small-scale rearrangements from the weaker scaffold set.
    """
    logger.info("Merging two iterations")

    #synchronizing scaffolds to the same permutations
    big_updated = _update_scaffolds(big_scaffolds, perm_container)
    small_updated = _update_scaffolds(small_scaffolds, perm_container)

    if rearrange:
        projector = RearrangementProjector(big_updated, small_updated, True)
        new_adj = projector.project()
        big_rearranged = build_scaffolds(new_adj, perm_container, False, False)
    else:
        big_rearranged = big_updated

    merged_scf = _merge_scaffolds(big_rearranged, small_updated)
    merged_scf = _merge_consecutive_contigs(merged_scf)

    if debugger.debugging:
        links_out = os.path.join(debugger.debug_dir, "merged.links")
        output_links(merged_scf, links_out)
        perms_out = os.path.join(debugger.debug_dir, "merged_scaffolds.txt")
        output_scaffolds_premutations(merged_scf, perms_out)

    return merged_scf


def get_breakpoints(scaffolds, bp_graph, perm_container):
    """
    Counts target-specific adjacencies in scaffolds
    """
    _updated_scaffolds = _update_scaffolds(scaffolds, perm_container)
    specific = 0
    for scf in scaffolds:
        for cnt in scf.contigs:
            for block_1, block_2 in cnt.perm.iter_pairs():
                genomes = bp_graph.genomes_support(-block_1.signed_id(),
                                                   block_2.signed_id())
                if set(genomes) == set([bp_graph.target]):
                    specific += 1

    logger.debug("Target-specific adjacencies in scaffolds: %d", specific)
    return specific


def _merge_consecutive_contigs(scaffolds):
    """
    Merges consecutive contig fragments originating from a same contig
    """
    new_scaffolds = []
    num_contigs = 0
    for scf in scaffolds:
        new_contigs = []

        cur_sign, cur_perm, cur_link = None, None, None
        for cnt in scf.contigs:
            consistent = False
            if cur_sign == cnt.sign and cnt.perm.chr_name == cur_perm.chr_name:
                if cur_sign > 0 and cur_perm.seq_end == cnt.perm.seq_start:
                    cur_perm.seq_end = cnt.perm.seq_end
                    cur_perm.blocks.extend(cnt.perm.blocks)
                    consistent = True
                if cur_sign < 0 and cur_perm.seq_start == cnt.perm.seq_end:
                    cur_perm.seq_start = cnt.perm.seq_start
                    cur_perm.blocks = cnt.perm.blocks + cur_perm.blocks
                    consistent = True

            if not consistent:
                if cur_perm:
                    new_contigs.append(Contig.with_perm(cur_perm, cur_sign, cur_link))
                cur_perm = deepcopy(cnt.perm)

            cur_sign = cnt.sign
            cur_link = cnt.link

        if cur_perm:
            new_contigs.append(Contig.with_perm(cur_perm, cur_sign, cur_link))
        num_contigs += len(new_contigs)
        new_scaffolds.append(Scaffold.with_contigs(scf.name, None,
                                                   None, new_contigs))

    logger.debug("Merging consequtive contigs: %d left", num_contigs)
    return new_scaffolds


def _update_scaffolds(scaffolds, perm_container):
    """
    Updates scaffolds wrt to given permutations
    """
    perm_index = defaultdict(list)
    for perm in perm_container.target_perms:
        perm_index[(perm.chr_name, perm.repeat_id)].append(perm)

    new_scaffolds = []
    for scf in scaffolds:
        new_contigs = []
        for contig in scf.contigs:
            inner_perms = []
            for new_perm in perm_index[(contig.perm.chr_name,
                                        contig.perm.repeat_id)]:
                if (contig.perm.seq_start <= new_perm.seq_start
                    < contig.perm.seq_end):
                    inner_perms.append(new_perm)
                    assert (contig.perm.seq_start < new_perm.seq_end
                            <= contig.perm.seq_end)

            if not inner_perms:
                logger.debug("Lost: %s", str(contig.perm))
                continue

            inner_perms.sort(key=lambda p: p.seq_start, reverse=contig.sign < 0)
            for prev_perm, next_perm in zip(inner_perms[:-1], inner_perms[1:]):
                if contig.sign > 0:
                    gap_length = next_perm.seq_start - prev_perm.seq_end
                else:
                    gap_length = prev_perm.seq_start - next_perm.seq_end
                support = [GenChrPair(prev_perm.genome_name, prev_perm.chr_name)]
                new_contigs.append(Contig.with_perm(prev_perm, contig.sign,
                                                    Link(gap_length, support)))
            new_contigs.append(Contig.with_perm(inner_perms[-1], contig.sign,
                                                copy(contig.link)))

        if len(new_contigs):
            new_scaffolds.append(Scaffold.with_contigs(scf.name, None,
                                                       None, new_contigs))
    return new_scaffolds


class RearrangementProjector:
    """
    This class handles the projection of rearrangements from weaker set
    of scaffolds and ensures that these rearrangements are small-scale
    """
    def __init__(self, old_scaffolds, new_scaffolds, conservative):
        self.old_scaffolds = old_scaffolds
        self.new_scaffolds = new_scaffolds
        self._build_bp_graph()
        self._build_adj_graph()
        self.conservative = conservative

    def connected_component_subgraphs(self,G):
        for c in nx.connected_components(G):
            yield G.subgraph(c)

    def project(self):
        #look for valid k-breaks
        num_kbreaks = 0
        #subgraphs = list(nx.connected_component_subgraphs(self.bp_graph))
        subgraphs = list(self.connected_component_subgraphs(self.bp_graph))
        for subgr in subgraphs:
            #this is a cycle
            if any(len(subgr[node]) != 2 for node in subgr.nodes):
                continue

            red_edges = []
            black_edges = []
            for (u, v, data) in subgr.edges(data=True):
                if data["scf_set"] == "old":
                    red_edges.append((u, v))
                else:
                    black_edges.append((u, v))

            if not self._good_k_break(red_edges, black_edges):
                continue

            num_kbreaks += 1
            for u, v in red_edges:
                self.bp_graph.remove_edge(u, v)
                self.adj_graph.remove_edge(u, v)
            for u, v in black_edges:
                #print(self.bp_graph[u][v])
                link = self.bp_graph[u][v][0]["link"]
                infinity = self.bp_graph[u][v][0]["infinity"]
                self.bp_graph.add_edge(u, v, scf_set="old",
                                       link=link, infinity=infinity)
                self.adj_graph.add_edge(u, v)

        logger.debug("Made %d k-breaks", num_kbreaks)
        adjacencies = {}
        for (u, v, data) in self.bp_graph.edges(data=True):
            if data["scf_set"] == "old":
                gap, support = 0, []
                if not data["infinity"]:
                    gap = data["link"].gap
                    support = data["link"].supporting_genomes
                adjacencies[u] = Adjacency(v, gap, support, data["infinity"])
                adjacencies[v] = Adjacency(u, gap, support, data["infinity"])

        return adjacencies

    def _good_k_break(self, old_edges, new_edges):
        """
        Checks that the break does not change chromomsome structure significantly
        """
        MIN_OVLP_SCORE = 0.9
        MAX_K_BREAK = 4
        if len(old_edges) > MAX_K_BREAK:
            return False

        new_adj_graph = self.adj_graph.copy()
        for u, v in old_edges:
            new_adj_graph.remove_edge(u, v)
        for u, v in new_edges:
            new_adj_graph.add_edge(u, v)

        old_sets = [set(g.nodes) for g in
                    #nx.connected_component_subgraphs(self.adj_graph)]
                    self.connected_component_subgraphs(self.adj_graph)]
        new_sets = [set(g.nodes) for g in
                    #nx.connected_component_subgraphs(new_adj_graph)]
                    self.connected_component_subgraphs(new_adj_graph)]
        if len(old_sets) != len(new_sets):
            return False

        for old_set in old_sets:
            max_overlap = 0
            best_score = 0
            for new_set in new_sets:
                overlap = len(old_set & new_set)
                if overlap > max_overlap:
                    max_overlap = overlap
                    best_score = float(overlap) / len(old_set | new_set)
            if best_score < MIN_OVLP_SCORE:
                return False

        return True

    def _build_bp_graph(self):
        """
        No repeats assumed!
        """
        old_contigs = set()
        for scf in self.old_scaffolds:
            for cnt in scf.contigs:
                old_contigs.add(cnt.name())

        ###creating 2-colored breakpoint graph
        bp_graph = nx.MultiGraph()
        for scf in self.old_scaffolds:
            for cnt_1, cnt_2 in zip(scf.contigs[:-1], scf.contigs[1:]):
                bp_graph.add_edge(cnt_1.right_end(), cnt_2.left_end(),
                                  scf_set="old", link=copy(cnt_1.link),
                                  scf_name=scf.name, infinity=False)
            #chromosome ends
            bp_graph.add_edge(scf.contigs[-1].right_end(),
                              scf.contigs[0].left_end(), scf_set="old",
                              infinity=True)

        for scf in self.new_scaffolds:
            prev_cont = None
            first_ctg = None
            pos = 0
            for pos, contig in enumerate(scf.contigs):
                if contig.name() in old_contigs:
                    prev_cont = deepcopy(contig)
                    first_ctg = prev_cont
                    break
            if prev_cont is None:
                continue

            for next_cont in scf.contigs[pos + 1:]:
                if next_cont.name() not in old_contigs:
                    prev_cont.link.gap += next_cont.length() + next_cont.link.gap
                    common_genomes = (set(prev_cont.link.supporting_genomes) &
                                      set(next_cont.link.supporting_genomes))
                    prev_cont.link.supporting_genomes = list(common_genomes)

                else:
                    bp_graph.add_edge(prev_cont.right_end(), next_cont.left_end(),
                                      scf_set="new", link=copy(prev_cont.link),
                                      scf_name=scf.name, infinity=False)
                    prev_cont = deepcopy(next_cont)

            bp_graph.add_edge(prev_cont.right_end(),
                              first_ctg.left_end(), scf_set="new",
                              infinity=True, link=None)

        self.bp_graph = bp_graph

    def _build_adj_graph(self):
        adj_graph = nx.Graph()
        for scf in self.old_scaffolds:
            for cnt_1, cnt_2 in zip(scf.contigs[:-1], scf.contigs[1:]):
                adj_graph.add_edge(cnt_1.right_end(), cnt_2.left_end())
            for cnt in scf.contigs:
                adj_graph.add_edge(cnt.left_end(), cnt.right_end())
            #chromosome ends
            adj_graph.add_edge(scf.contigs[-1].right_end(),
                               scf.contigs[0].left_end())
        self.adj_graph = adj_graph


def _merge_scaffolds(big_scaffolds, small_scaffolds):
    """
    Performs the final merging step
    """
    count_diff_scaf = 0
    count_diff_orient = 0
    count_inconsistent = 0

    total_success = 0
    total_fail = 0
    total_inserted = 0
    not_found = 0

    big_count = defaultdict(int)
    for scf in big_scaffolds:
        for c in scf.contigs:
            big_count[c.perm] += 1

    small_count = defaultdict(int)
    for scf in small_scaffolds:
        for c in scf.contigs:
            small_count[c.perm] += 1

    repeats = set(seq for (seq, count) in
                  chain(list(big_count.items()), list(small_count.items())) if count > 1)
    big_unique = set(seq for (seq, count) in big_count.items() if count == 1)

    small_index = {}
    for scf in small_scaffolds:
        for pos, contig in enumerate(scf.contigs):
            if contig.perm not in repeats:
                assert contig.perm not in small_index
                small_index[contig.perm] = (scf, pos)

    new_scafflods = []
    for big_scf in big_scaffolds:
        new_contigs = []
        #non_repeats = list(filter(lambda i: big_scf.contigs[i].perm
        #                                not in repeats,
        #                          xrange(len(big_scf.contigs))))
        non_repeats = [i for i in range(len(big_scf.contigs))
                       if big_scf.contigs[i].perm not in repeats]
        for left_idx, right_idx in zip(non_repeats[:-1], non_repeats[1:]):
            left_cnt = big_scf.contigs[left_idx]
            right_cnt = big_scf.contigs[right_idx]

            consistent = False
            weak_contigs = None
            link_to_change = None
            if (left_cnt.perm in small_index and
                right_cnt.perm in small_index):
                consistent = True
                left_scf, left_pos = small_index[left_cnt.perm]
                right_scf, right_pos = small_index[right_cnt.perm]

                big_sign = left_cnt.sign == right_cnt.sign
                small_sign = (left_scf.contigs[left_pos].sign ==
                              right_scf.contigs[right_pos].sign)

                if left_scf != right_scf:
                    count_diff_scaf += 1
                    consistent = False
                elif big_sign != small_sign:
                    count_diff_orient += 1
                    consistent = False
                else:
                    same_dir = left_pos < right_pos
                    if not same_dir:
                        left_pos, right_pos = right_pos, left_pos

                    weak_contigs = left_scf.contigs[left_pos + 1 : right_pos]
                    if any(c.perm in big_unique for c in weak_contigs):
                        count_inconsistent += 1
                        consistent = False

                    link_to_change = copy(left_scf.contigs[left_pos].link)
                    #reverse complement
                    if weak_contigs and not same_dir:
                        link_to_change = copy(left_scf.contigs[right_pos - 1].link)
                        weak_contigs = [c.reverse_copy() for c in weak_contigs[::-1]]
                        for pw, nw in zip(weak_contigs[:-1], weak_contigs[1:]):
                            pw.link = copy(nw.link)
                        weak_contigs[-1].link = copy(left_scf.contigs[left_pos].link)

            else:
                not_found += 1

            new_contigs.append(left_cnt)
            if consistent and weak_contigs:
                new_contigs[-1].link = link_to_change
                new_contigs.extend(weak_contigs)
                total_success += 1
                total_inserted += len(weak_contigs)
                #logger.debug("Inserting '{0}' between {1} and {2}"
                #             .format(map(lambda c: c.perm, weak_contigs),
                #                     left_cnt, right_cnt))
            else:
                new_contigs.extend(big_scf.contigs[left_idx+1:right_idx])
                total_fail += 1

        if len(new_contigs) > 1:
            new_contigs.append(right_cnt)
            s = Scaffold(big_scf.name)
            s.contigs = new_contigs
            new_scafflods.append(s)
        else:   #because of repeats
            new_scafflods.append(big_scf)

    logger.debug("Fail: not found: %d", not_found)
    logger.debug("Fail: different scaffolds: %d", count_diff_scaf)
    logger.debug("Fail: different orientatilns: %d", count_diff_orient)
    logger.debug("Fail: inconsistent: %d", count_inconsistent)
    logger.debug("Total success: %d", total_success)
    logger.debug("Total fail: %d", total_fail)
    logger.debug("Total inserted: %d", total_inserted)

    num_contigs = 0
    for scf in new_scafflods:
        num_contigs += len(scf.contigs)
    logger.debug("Result: %d contigs in %d scaffolds", num_contigs, len(new_scafflods))

    return new_scafflods
