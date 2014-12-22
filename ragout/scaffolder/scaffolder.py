#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module assembles contigs into scaffolds with respect
to given adjacencies. Also, it outputs scaffolds in different
formats
"""

from collections import defaultdict
from itertools import repeat
import os
import copy
import logging

from ragout.shared.debug import DebugConfig
from ragout.shared.datatypes import Contig, Scaffold, Link
from .output_generator import output_links

logger = logging.getLogger()
debugger = DebugConfig.get_instance()


def get_scaffolds(adjacencies, perm_container):
    """
    Assembles scaffolds
    """
    logger.info("Building scaffolds")
    contigs, contig_index = _make_contigs(perm_container)
    scaffolds = _extend_scaffolds(adjacencies, contigs, contig_index)
    scaffolds = list(filter(lambda s: len(s.contigs) > 1, scaffolds))

    if debugger.debugging:
        links_out = os.path.join(debugger.debug_dir, "scaffolds.links")
        output_links(scaffolds, links_out)
        perms_out = os.path.join(debugger.debug_dir, "scaffolds_perms.txt")
        _output_scaffold_premutations(scaffolds, perms_out)

    return scaffolds


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
        scf = Scaffold.with_contigs(scf_name, contig.blocks[0],
                                    contig.blocks[-1], [contig])
        scaffolds.append(scf)

        #go right
        while scf.right in adjacencies:
            adj_block = adjacencies[scf.right].block
            adj_distance = adjacencies[scf.right].distance
            adj_supporting_genomes = adjacencies[scf.right].supporting_genomes
            assert len(contig_index[abs(adj_block)]) == 1

            contig = contig_index[abs(adj_block)][0]
            if contig in visited:
                break

            if adj_block in [contig.blocks[0], -contig.blocks[-1]]:
                scf.contigs[-1].link = Link(adj_distance, adj_supporting_genomes)
                scf.contigs.append(contig)
                visited.add(contig)

                if contig.blocks[0] == adj_block:
                    scf.right = contig.blocks[-1]
                else:
                    scf.contigs[-1].sign = -1
                    scf.right = -contig.blocks[0]

                continue

            break

        #go left
        while -scf.left in adjacencies:
            adj_block = -adjacencies[-scf.left].block
            adj_distance = adjacencies[-scf.left].distance
            adj_supporting_genomes = adjacencies[-scf.left].supporting_genomes
            assert len(contig_index[abs(adj_block)]) == 1

            contig = contig_index[abs(adj_block)][0]
            if contig in visited:
                break

            if adj_block in [contig.blocks[-1], -contig.blocks[0]]:
                scf.contigs.insert(0, contig)
                scf.contigs[0].link = Link(adj_distance, adj_supporting_genomes)
                visited.add(contig)

                if contig.blocks[-1] == adj_block:
                    scf.left = contig.blocks[0]

                else:
                    scf.contigs[0].sign = -1
                    scf.left = -contig.blocks[-1]

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
                if contig.sign > 0:
                    blocks.extend(contig.blocks)
                else:
                    rev_compl = map(lambda b: -b, contig.blocks[::-1])
                    blocks.extend(rev_compl)

            f.write(">" + scf.name + "\n")
            for block in blocks:
                f.write("{0:+} ".format(block))
            f.write("$\n")


def _make_contigs(perm_container):
    """
    Converts permutations into contigs
    """
    contigs = []
    index = defaultdict(list)
    for perm in perm_container.target_perms:
        assert len(perm.blocks)

        contigs.append(Contig(perm.chr_name))
        for block in perm.blocks:
            index[block.block_id].append(contigs[-1])
            contigs[-1].blocks.append(block.signed_id())

    return contigs, index
