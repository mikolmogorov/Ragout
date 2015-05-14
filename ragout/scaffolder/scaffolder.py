#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module assembles contigs into scaffolds with respect
to given adjacencies. Also, it outputs scaffolds in different
formats
"""

from collections import defaultdict, namedtuple
from itertools import repeat
import os
import copy
import logging

import networkx as nx

from ragout.shared.debug import DebugConfig
from ragout.shared.datatypes import Contig, Scaffold, Link
from .output_generator import output_links

logger = logging.getLogger()
debugger = DebugConfig.get_instance()

Adjacency = namedtuple("Adjacency", ["block", "distance", "supporting_genomes"])

def get_scaffolds(adjacencies, perm_container):
    """
    Assembles scaffolds
    """
    logger.info("Building scaffolds")
    contigs, contig_index = _make_contigs(perm_container)
    scaffolds = _extend_scaffolds(adjacencies, contigs, contig_index)
    scaffolds = list(filter(lambda s: len(s.contigs) > 1, scaffolds))

    num_contigs = sum(map(lambda s: len(s.contigs), scaffolds))
    logger.debug("{0} contigs were joined into {1} scaffolds"
                        .format(num_contigs, len(scaffolds)))

    if debugger.debugging:
        links_out = os.path.join(debugger.debug_dir, "scaffolds.links")
        output_links(scaffolds, links_out)
        perms_out = os.path.join(debugger.debug_dir, "scaffolds_perms.txt")
        _output_scaffold_premutations(scaffolds, perms_out)

    return scaffolds


def update_scaffolds(scaffolds, new_perm_container):
    """
    Updates scaffolds wrt to new permutations
    """
    by_chr_name = defaultdict(list)
    for perm in new_perm_container.target_perms:
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


def project_rearrangements(old_scaffolds, new_scaffolds):
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


def _extend_scaffolds(adjacencies, contigs, contig_index):
    """
    Assembles contigs into scaffolds
    """
    scaffolds = []
    visited = set()
    counter = [0]

    def extend_scaffold(contig):
        visited.add(contig)
        scf_name = "ragout-scaffold-{0}".format(counter[0])
        counter[0] += 1
        scf = Scaffold.with_contigs(scf_name, contig.left_end(),
                                    contig.right_end(), [contig])
        scaffolds.append(scf)

        #go right
        while scf.right in adjacencies:
            adj_block = adjacencies[scf.right].block
            adj_distance = adjacencies[scf.right].distance
            adj_supporting_genomes = adjacencies[scf.right].supporting_genomes

            contig = contig_index[abs(adj_block)]
            if contig in visited:
                break

            if adj_block in [contig.left_end(), contig.right_end()]:
                if contig.left_end() == adj_block:
                    scf.contigs.append(contig)
                else:
                    scf.contigs.append(contig.reverse_copy())

                flank = scf.contigs[-2].right_gap() + scf.contigs[-1].left_gap()
                gap = max(0, adj_distance - flank)
                scf.contigs[-2].link = Link(gap, adj_supporting_genomes)

                scf.right = scf.contigs[-1].right_end()
                visited.add(contig)
                continue

            break

        #go left
        while scf.left in adjacencies:
            adj_block = adjacencies[scf.left].block
            adj_distance = adjacencies[scf.left].distance
            adj_supporting_genomes = adjacencies[scf.left].supporting_genomes

            contig = contig_index[abs(adj_block)]
            if contig in visited:
                break

            if adj_block in [contig.right_end(), contig.left_end()]:
                if contig.right_end() == adj_block:
                    scf.contigs.insert(0, contig)
                else:
                    scf.contigs.insert(0, contig.reverse_copy())

                flank = scf.contigs[0].right_gap() + scf.contigs[1].left_gap()
                gap = max(0, adj_distance - flank)
                scf.contigs[0].link = Link(gap, adj_supporting_genomes)

                scf.left = scf.contigs[0].left_end()
                visited.add(contig)
                continue

            break

    for contig in contigs:
        if contig not in visited:
            extend_scaffold(contig)

    return scaffolds


def _output_scaffold_premutations(scaffolds, out_file):
    with open(out_file, "w") as f:
        permutations = []
        for scf in scaffolds:
            blocks = []
            for contig in scf.contigs:
                blocks.extend(contig.signed_perm())

            f.write(">" + scf.name + "\n")
            for block in blocks:
                f.write("{0:+} ".format(block))
            f.write("$\n")


def _make_contigs(perm_container):
    """
    Converts permutations into contigs
    """
    contigs = []
    index = {}
    for perm in perm_container.target_perms:
        assert len(perm.blocks)

        contigs.append(Contig(perm.name(), perm))
        for block in perm.blocks:
            assert block.block_id not in index
            index[block.block_id] = contigs[-1]

    return contigs, index
