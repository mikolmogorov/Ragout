#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module assembles contigs into scaffolds with respect
to given adjacencies. Also, it outputs scaffolds in different
formats
"""

from __future__ import absolute_import
from __future__ import division
from collections import defaultdict
import os
import logging

from ragout.shared.debug import DebugConfig
from ragout.shared.datatypes import (Contig, Scaffold, Link,
                                     output_scaffolds_premutations,
                                     output_permutations)
from ragout.scaffolder.output_generator import output_links
from ragout.breakpoint_graph.inferer import Adjacency
from ragout.six.moves import range, zip


logger = logging.getLogger()
debugger = DebugConfig.get_instance()


def build_scaffolds(adjacencies, perm_container, debug_output,
                    correct_distances):
    """
    Assembles scaffolds wrt to inferred adjacencies
    """
    if debug_output:
        logger.debug("Building scaffolds")
    contigs, contig_index = _make_contigs(perm_container)
    scaffolds = _extend_scaffolds(adjacencies, contigs, contig_index,
                                  correct_distances)
    num_contigs = sum([len(s.contigs) for s in scaffolds])
    logger.debug("%d contigs were joined into %d scaffolds",
                 num_contigs, len(scaffolds))

    if debugger.debugging and debug_output:
        links_out = os.path.join(debugger.debug_dir, "scaffolder.links")
        output_links(scaffolds, links_out)
        contigs_out = os.path.join(debugger.debug_dir, "scaffolder_contigs.txt")
        output_permutations(perm_container.target_perms, contigs_out)
        perms_out = os.path.join(debugger.debug_dir, "scaffolder_scaffolds.txt")
        output_scaffolds_premutations(scaffolds, perms_out)

    return scaffolds


def update_gaps(scaffolds):
    """
    Do it in the very end
    """
    for scf in scaffolds:
        for c1, c2 in zip(scf.contigs[:-1], scf.contigs[1:]):
            c1.link.gap -= c1.right_gap() + c2.left_gap()


def assign_scaffold_names(scaffolds, perm_container, ref_genome):
    """
    Names scaffolds according to homology to a chosen reference genome.
    Also ensures that scaffolds and corresponding reference chromosomes
    have the same strand.
    """
    MIN_RATE = 0.1
    PREFIX = "chr"
    chr_index = {}
    for perm in perm_container.ref_perms:
        if perm.genome_name == ref_genome:
            for block in perm.blocks:
                chr_index[block.block_id] = perm.chr_name, block.sign

    assigned_names = {}
    need_rev_compl = {}
    for scf in scaffolds:
        scf_index = defaultdict(int)
        sign_agreement = 0
        total = 0
        for contig in scf.contigs:
            for block in contig.perm.blocks:
                if block.block_id in chr_index:
                    chrom, sign = chr_index[block.block_id]
                    scf_index[chrom] += 1
                    total += 1
                    sign_agreement += int(sign == block.sign * contig.sign)

        name_str = PREFIX
        for chrom in sorted(scf_index, key=scf_index.get, reverse=True):
            if scf_index[chrom] > MIN_RATE * total:
                name_str += "_" + chrom
            else:
                break
        assigned_names[scf] = name_str
        need_rev_compl[scf] = sign_agreement < total // 2

    #in case of same names
    same_names = defaultdict(list)
    for scf, name in assigned_names.items():
        same_names[name].append(scf)
    for name, scf_list in same_names.items():
        scf_list.sort(key=lambda s: len(s.contigs), reverse=True)
        unlocalized = scf_list[1:]
        for scf in unlocalized:
            assigned_names[scf] += "_unlocalized"
        if len(unlocalized) > 1:
            for num, scf in enumerate(unlocalized):
                assigned_names[scf] += "." + str(num + 1)

    for scf in scaffolds:
        scf.name = assigned_names[scf]
        if need_rev_compl[scf]:
            new_contigs = [c.reverse_copy() for c in scf.contigs][::-1]
            for i in range(len(new_contigs) - 1):
                new_contigs[i].link = new_contigs[i + 1].link
            new_contigs[-1].link = Link(0, [])
            scf.contigs = new_contigs


def _extend_scaffolds(adjacencies, contigs, contig_index, correct_distances):
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

        already_complete = (scf.right in adjacencies and
                            adjacencies[scf.right].block == scf.left and
                            adjacencies[scf.right].infinity)
        if already_complete:
            scaffolds.append(scf)
            return

        #go right
        while scf.right in adjacencies and not adjacencies[scf.right].infinity:
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
                gap = adj_distance - flank if correct_distances else adj_distance
                scf.contigs[-2].link = Link(gap, adj_supporting_genomes)

                scf.right = scf.contigs[-1].right_end()
                visited.add(contig)
                continue

            break

        #go left
        while scf.left in adjacencies and not adjacencies[scf.left].infinity:
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
                gap = adj_distance - flank if correct_distances else adj_distance
                scf.contigs[0].link = Link(gap, adj_supporting_genomes)

                scf.left = scf.contigs[0].left_end()
                visited.add(contig)
                continue

            break

        if len(scf.contigs) > 1:
            scaffolds.append(scf)

    #for contig in contigs:
    for contig in sorted(contigs, key=lambda c: c.signed_name()):
        if contig not in visited:
            extend_scaffold(contig)

    return scaffolds


def _make_contigs(perm_container):
    """
    A helper function to make Contig structures
    """
    contigs = []
    index = {}
    for perm in perm_container.target_perms:
        assert len(perm.blocks)
        contigs.append(Contig.with_perm(perm))
        for block in perm.blocks:
            assert block.block_id not in index
            index[block.block_id] = contigs[-1]

    return contigs, index
