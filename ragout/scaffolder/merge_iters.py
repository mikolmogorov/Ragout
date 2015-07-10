#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module provides some functions for
moving between two consecutive iterations
"""

from collections import namedtuple, defaultdict
from itertools import product, chain
import os
import logging
from copy import deepcopy

import networkx as nx

from ragout.shared.debug import DebugConfig
from ragout.shared.datatypes import (ContigWithPerm, Scaffold, Permutation, Link,
                                     output_scaffolds_premutations, output_permutations)
from ragout.scaffolder.output_generator import output_links
from ragout.scaffolder.scaffolder import build_scaffolds


logger = logging.getLogger()
debugger = DebugConfig.get_instance()
Adjacency = namedtuple("Adjacency", ["block", "distance", "supporting_genomes"])


def merge_scaffolds(big_scaffolds, small_scaffolds, perm_container, rearrange):
    logger.info("Merging two iterations")
    #synchronizing scaffolds
    big_updated = _update_scaffolds(big_scaffolds, perm_container)
    small_updated = _update_scaffolds(small_scaffolds, perm_container)

    if rearrange:
        new_adj = _project_rearrangements(big_updated, small_updated)
        big_rearranged = build_scaffolds(new_adj, perm_container, False, False)
    else:
        big_rearranged = big_updated

    merged_scf = _merge_scaffolds(big_rearranged, small_updated)
    merged_scf = _merge_consecutive_contigs(merged_scf)

    if debugger.debugging:
        links_out = os.path.join(debugger.debug_dir, "merged.links")
        output_links(merged_scf, links_out)
        #contigs_out = os.path.join(debugger.debug_dir, "merged_contigs.txt")
        #output_permutations(perm_container.target_perms, contigs_out)
        perms_out = os.path.join(debugger.debug_dir, "merged_scaffolds.txt")
        output_scaffolds_premutations(merged_scf, perms_out)

    return merged_scf


def refine_scaffolds(scaffolds, adj_refiner, perm_container):
    updated_scf = _update_scaffolds(scaffolds, perm_container)
    adj = adj_refiner.refine_adjacencies(updated_scf)
    new_scf = build_scaffolds(adj, perm_container)
    return _merge_consecutive_contigs(new_scf)


def _merge_consecutive_contigs(scaffolds):
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
                    new_contigs.append(ContigWithPerm(cur_perm, cur_sign, cur_link))
                cur_perm = deepcopy(cnt.perm)

            cur_sign = cnt.sign
            cur_link = cnt.link

        if cur_perm:
            new_contigs.append(ContigWithPerm(cur_perm, cur_sign, cur_link))
        num_contigs += len(new_contigs)
        new_scaffolds.append(Scaffold.with_contigs(scf.name, None,
                                                   None, new_contigs))

    logger.debug("Merging consequtive contigs: {0} left".format(num_contigs))
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
                logger.debug("Lost: {0}".format(contig.perm))
            inner_perms.sort(key=lambda p: p.seq_start, reverse=contig.sign < 0)
            for new_perm in inner_perms:
                new_contigs.append(ContigWithPerm(new_perm, contig.sign,
                                                  Link(0, ["^_^"])))
            new_contigs[-1].link = contig.link

        new_scaffolds.append(Scaffold.with_contigs(scf.name, None,
                                                   None, new_contigs))
    return new_scaffolds


def _project_rearrangements(old_scaffolds, new_scaffolds):
    """
    No repeats assumed!
    """
    old_contigs = set()
    for scf in old_scaffolds:
        for cnt in scf.contigs:
            old_contigs.add(cnt.name())

    ###creating 2-colored breakpoint graph
    bp_graph = nx.MultiGraph()
    for scf in old_scaffolds:
        for cnt_1, cnt_2 in zip(scf.contigs[:-1], scf.contigs[1:]):
            bp_graph.add_edge(cnt_1.right_end(), cnt_2.left_end(),
                              scf_set="old", link=cnt_1.link,
                              scf_name=scf.name, c1=cnt_1, c2=cnt_2)

    for scf in new_scaffolds:
        prev_cont = None
        for pos, contig in enumerate(scf.contigs):
            if contig.name() in old_contigs:
                prev_cont = contig
                break
        if prev_cont is None:
            continue

        for next_cont in scf.contigs[pos + 1:]:
            if next_cont.name() not in old_contigs:
                prev_cont.link.gap += next_cont.length() + next_cont.link.gap
                common_genomes = (set(prev_cont.link.supporting_genomes) &
                                  set(next_cont.link.supporting_genomes))
                prev_cont.link.supporting_genomes = list(common_genomes)
                continue

            bp_graph.add_edge(prev_cont.right_end(), next_cont.left_end(),
                              scf_set="new", link=prev_cont.link,
                              c1=prev_cont, c2=next_cont)
            prev_cont = next_cont
    ###

    #now look for valid k-breaks
    num_kbreaks = 0
    subgraphs = list(nx.connected_component_subgraphs(bp_graph))
    for subgr in subgraphs:
        #this is a cycle
        if any(len(subgr.neighbors(node)) != 2 for node in subgr.nodes()):
            continue

        red_edges = []
        black_edges = []
        scaffolds_involved = set()
        #logger.debug(">>>k-break")

        num_red = 0
        for (u, v, data) in subgr.edges_iter(data=True):
            if data["scf_set"] == "old":
                red_edges.append((u, v))
                scaffolds_involved.add(data["scf_name"])
            else:
                black_edges.append((u, v))
                #logger.debug("{0} -- {1}".format(data["c1"].signed_name(),
                #                                 data["c2"].signed_name()))
                num_red += int(_red_supported(data["c1"], data["c2"]))

        if num_red != len(subgr) / 2 - 1:
            continue

        assert len(red_edges) == len(black_edges)

        #if len(subgr) > 3:
        #    continue
        #if len(scaffolds_involved) > 2:
        #    continue
        #logger.debug("{0}-break in {1} scaffolds".format(len(subgr) / 2,
        #                                                 len(scaffolds_involved)))
        num_kbreaks += 1

        for u, v in red_edges:
            bp_graph.remove_edge(u, v)
        for u, v in black_edges:
            link = bp_graph[u][v][0]["link"]
            bp_graph.add_edge(u, v, scf_set="old", link=link)

    logger.debug("Made {0} k-breaks".format(num_kbreaks))
    adjacencies = {}
    for (u, v, data) in bp_graph.edges_iter(data=True):
        if data["scf_set"] == "old":
            adjacencies[u] = Adjacency(v, data["link"].gap,
                                       data["link"].supporting_genomes)
            adjacencies[v] = Adjacency(u, data["link"].gap,
                                       data["link"].supporting_genomes)

    return adjacencies


def _red_supported(ctg_1, ctg_2):
    if ctg_1.perm.chr_name != ctg_2.perm.chr_name:
        return False

    left = ctg_1.perm.seq_end if ctg_1.sign > 0 else ctg_1.perm.seq_start
    right = ctg_2.perm.seq_start if ctg_2.sign > 0 else ctg_2.perm.seq_end
    return left == right


def _merge_scaffolds(big_scaffolds, small_scaffolds):
    """
    Merges two assemblies assuming that big one is correct
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
                  chain(big_count.items(), small_count.items()) if count > 1)
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
        non_repeats = list(filter(lambda i: big_scf.contigs[i].perm
                                        not in repeats,
                                  xrange(len(big_scf.contigs))))
        for left_idx, right_idx in zip(non_repeats[:-1], non_repeats[1:]):
            left_cnt = big_scf.contigs[left_idx]
            right_cnt = big_scf.contigs[right_idx]

            consistent = False
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

                    if not same_dir:
                        weak_contigs = list(map(lambda c: c.reverse_copy(),
                                                weak_contigs[::-1]))
                    link_to_change = left_scf.contigs[left_pos].link
            else:
                not_found += 1

            new_contigs.append(left_cnt)
            if consistent:
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

    logger.debug("Fail: not found: {0}".format(not_found))
    logger.debug("Fail: different scaffolds: {0}".format(count_diff_scaf))
    logger.debug("Fail: different orientatilns: {0}".format(count_diff_orient))
    logger.debug("Fail: inconsistent: {0}".format(count_inconsistent))
    logger.debug("Total success: {0}".format(total_success))
    logger.debug("Total fail: {0}".format(total_fail))
    logger.debug("Total inserted: {0}".format(total_inserted))

    num_contigs = 0
    for scf in new_scafflods:
        num_contigs += len(scf.contigs)
    logger.debug("Result: {0} contigs in {1} scaffolds"
                                    .format(num_contigs, len(new_scafflods)))

    return new_scafflods
