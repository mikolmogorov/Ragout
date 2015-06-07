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

from ragout.shared.debug import DebugConfig
from ragout.shared.datatypes import (ContigWithPerm, Scaffold, Link,
                                     output_scaffolds_premutations,
                                     output_permutations)
from ragout.scaffolder.output_generator import output_links


logger = logging.getLogger()
debugger = DebugConfig.get_instance()


def build_scaffolds(adjacencies, perm_container, debug_output=True,
                    correct_distances=True):
    """
    Assembles scaffolds
    """
    logger.info("Building scaffolds")
    contigs, contig_index = _make_contigs(perm_container)
    scaffolds = _extend_scaffolds(adjacencies, contigs, contig_index,
                                  correct_distances)
    scaffolds = list(filter(lambda s: len(s.contigs) > 1, scaffolds))

    num_contigs = sum(map(lambda s: len(s.contigs), scaffolds))
    logger.debug("{0} contigs were joined into {1} scaffolds"
                        .format(num_contigs, len(scaffolds)))

    if debugger.debugging and debug_output:
        links_out = os.path.join(debugger.debug_dir, "scaffolder.links")
        output_links(scaffolds, links_out)
        contigs_out = os.path.join(debugger.debug_dir, "scaffolder_contigs.txt")
        output_permutations(perm_container.target_perms, contigs_out)
        perms_out = os.path.join(debugger.debug_dir, "scaffolder_scaffolds.txt")
        output_scaffolds_premutations(scaffolds, perms_out)

    return scaffolds


def assign_scaffold_names(scaffolds, perm_container, ref_genome):
    MIN_RATE = 0.1
    PREFIX = "pseudochr"
    chr_index = {}
    for perm in perm_container.ref_perms:
        if perm.genome_name == ref_genome:
            for block in perm.blocks:
                chr_index[block.block_id] = perm.chr_name

    assigned_names = {}
    for scf in scaffolds:
        scf_index = defaultdict(int)
        total = 0
        for contig in scf.contigs:
            for block in contig.perm.blocks:
                if block.block_id in chr_index:
                    scf_index[chr_index[block.block_id]] += 1
                    total += 1

        name_str = PREFIX
        for chrom in sorted(scf_index, key=scf_index.get, reverse=True):
            if scf_index[chrom] > MIN_RATE * total:
                name_str += "." + chrom
            else:
                break
        assigned_names[scf] = name_str

    same_names = defaultdict(list)
    for scf, name in assigned_names.items():
        same_names[name].append(scf)
    for name, scf_list in same_names.items():
        if len(scf_list) == 1:
            continue
        for num, scf in enumerate(scf_list):
            assigned_names[scf] += "." + str(num + 1)

    for scf in scaffolds:
        scf.name = assigned_names[scf]


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
                gap = adj_distance - flank if correct_distances else adj_distance
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
                gap = adj_distance - flank if correct_distances else adj_distance
                scf.contigs[0].link = Link(gap, adj_supporting_genomes)

                scf.left = scf.contigs[0].left_end()
                visited.add(contig)
                continue

            break

    for contig in contigs:
        if contig not in visited:
            extend_scaffold(contig)

    return scaffolds


def _make_contigs(perm_container):
    """
    Converts permutations into contigs
    """
    contigs = []
    index = {}
    for perm in perm_container.target_perms:
        assert len(perm.blocks)
        contigs.append(ContigWithPerm(perm))
        for block in perm.blocks:
            assert block.block_id not in index
            index[block.block_id] = contigs[-1]

    return contigs, index
