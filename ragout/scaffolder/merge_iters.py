#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module provides some functions for
moving between two consecutive iterations
"""

from collections import namedtuple, defaultdict
from itertools import product, chain
import sys
import logging

import networkx as nx

from ragout.shared.datatypes import Contig, Scaffold
from ragout.scaffolder.scaffolder import build_scaffolds

logger = logging.getLogger()
Adjacency = namedtuple("Adjacency", ["block", "distance", "supporting_genomes"])


def merge_scaffolds(big_scaffolds, small_scaffolds, perm_container, rearrange):
    logger.info("Merging two iterations")
    big_updated = _update_scaffolds(big_scaffolds, perm_container)
    small_updated = _update_scaffolds(small_scaffolds, perm_container)

    if rearrange:
        #TODO: keep only common permutations
        new_adj = _project_rearrangements(big_updated, small_updated)
        big_rearranged = build_scaffolds(new_adj, perm_container)
    else:
        big_rearranged = big_updated

    return _merge(big_rearranged, small_updated)


def _update_scaffolds(scaffolds, perm_container):
    """
    Updates scaffolds wrt to given permutations
    """
    by_chr_name = defaultdict(list)
    for perm in perm_container.target_perms:
        by_chr_name[perm.chr_name].append(perm)

    new_scaffolds = []
    for scf in scaffolds:
        new_contigs = []
        for contig in scf.contigs:
            inner_perms = []
            for new_perm in by_chr_name[contig.perm.chr_name]:
                if (contig.perm.seq_start <= new_perm.seq_start
                    < contig.perm.seq_end):
                    inner_perms.append(new_perm)
                    assert (contig.perm.seq_start < new_perm.seq_end
                            <= contig.perm.seq_end)

            if not inner_perms:
                logger.debug("Lost: {0}".format(contig.perm))
            inner_perms.sort(key=lambda p: p.seq_start, reverse=contig.sign < 0)
            for new_perm in inner_perms:
                new_contigs.append(Contig(new_perm.name(), new_perm,
                                          contig.sign))

        new_scaffolds.append(Scaffold.with_contigs(scf.name, None,
                                                   None, new_contigs))
    return new_scaffolds


def _project_rearrangements(old_scaffolds, new_scaffolds):
    bp_graph = nx.MultiGraph()
    for scf in old_scaffolds:
        for cnt_1, cnt_2 in zip(scf.contigs[:-1], scf.contigs[1:]):
            bp_graph.add_edge(cnt_1.right_end(), cnt_2.left_end(),
                              scf_set="old", scf_name=scf.name)
    for scf in new_scaffolds:
        for cnt_1, cnt_2 in zip(scf.contigs[:-1], scf.contigs[1:]):
            bp_graph.add_edge(cnt_1.right_end(), cnt_2.left_end(),
                              scf_set="new", scf_name=scf.name)

    #now look for valid 2-breaks
    subgraphs = list(nx.connected_component_subgraphs(bp_graph))
    for subgr in subgraphs:
        if len(subgr) != 4:
            continue

        #this is a cycle
        if any(len(subgr.neighbors(node)) != 2 for node in subgr.nodes()):
            continue

        red_edges = []
        black_edges = []
        for (u, v, data) in subgr.edges_iter(data=True):
            if data["scf_set"] == "old":
                red_edges.append((u, v))
            else:
                black_edges.append((u, v))
        assert len(red_edges) == 2 and len(black_edges) == 2
        logger.debug("2-break!")

        for u, v in red_edges:
            bp_graph.remove_edge(u, v)
        for u, v in black_edges:
            bp_graph.add_edge(u, v, scf_set="old")

    adjacencies = {}
    for (u, v, data) in bp_graph.edges_iter(data=True):
        if data["scf_set"] == "old":
            adjacencies[u] = Adjacency(v, 0, [])
            adjacencies[v] = Adjacency(u, 0, [])

    return adjacencies


def _merge(big_scaffolds, small_scaffolds):
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
