#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module assembles contigs into scaffolds with respect
to given adjacencies. Also, it outputs scaffolds in different
formats
"""

from collections import defaultdict
import os
import copy
import logging

from ragout.shared.debug import DebugConfig
from ragout.shared.datatypes import Contig, Scaffold
from ragout.parsers.fasta_parser import write_fasta_dict, reverse_complement

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
        file_out = os.path.join(debugger.debug_dir, "scaffolds.ord")
        output_order(scaffolds, file_out)

    return scaffolds


def output_order(scaffolds, out_order):
    """
    Outputs scaffolds to file in "ord" format
    """
    out_order_stream = open(out_order, "w")
    for scf in scaffolds:
        out_order_stream.write(">" + scf.name + "\n")
        for contig in scf.contigs:
            out_order_stream.write(str(contig) + "\n")


def output_fasta(contigs_fasta, scaffolds, out_file):
    """
    Outputs scaffodls to file in "fasta" format
    """
    logger.info("Generating FASTA output")
    used_contigs = set()
    out_fasta_dict = {}

    scf_length = []
    for scf in scaffolds:
        scf_seqs = []
        for contig in scf.contigs:
            cont_seq = contigs_fasta[contig.name]
            used_contigs.add(contig.name)

            if contig.sign < 0:
                cont_seq = reverse_complement(cont_seq)

            if contig.gap >= 0:
                scf_seqs.append(cont_seq)
                scf_seqs.append("N" * contig.gap)
            else:
                scf_seqs.append(cont_seq[:contig.gap])

        scf_seq = "".join(scf_seqs)
        scf_length.append(len(scf_seq))
        out_fasta_dict[scf.name] = scf_seq
    write_fasta_dict(out_fasta_dict, out_file)

    #add some statistics
    used_count = 0
    used_len = 0
    unused_count = 0
    unused_len = 0
    for h in contigs_fasta:
        if h in used_contigs:
            used_count += 1
            used_len += len(contigs_fasta[h])
        else:
            unused_count += 1
            unused_len += len(contigs_fasta[h])
    assembly_len = unused_len + used_len
    used_perc = 100 * float(used_len) / assembly_len
    unused_perc = 100 * float(unused_len) / assembly_len
    contigs_length = [len(c) for c in contigs_fasta.values()]

    logger.info("Assembly statistics:\n\n"
                "\tScaffolds count:\t{0}\n"
                "\tUsed contigs count:\t{1}\n"
                "\tUsed contigs length:\t{2} ({3:2.4}%)\n"
                "\tUnused contigs count:\t{4}\n"
                "\tUnused contigs length:\t{5} ({6:2.4}%)\n"
                "\tContigs N50: \t\t{7}\n"
                "\tScaffolds N50:\t\t{8}\n"
                .format(len(scaffolds), used_count, used_len, used_perc,
                        unused_count, unused_len, unused_perc,
                        _calc_n50(contigs_length, unused_len + used_len),
                        _calc_n50(scf_length, unused_len + used_len)))


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
            assert len(contig_index[abs(adj_block)]) == 1

            contig = contig_index[abs(adj_block)][0]
            if contig in visited:
                break

            if adj_block in [contig.blocks[0], -contig.blocks[-1]]:
                scf.contigs[-1].gap = adj_distance
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
            assert len(contig_index[abs(adj_block)]) == 1

            contig = contig_index[abs(adj_block)][0]
            if contig in visited:
                break

            if adj_block in [contig.blocks[-1], -contig.blocks[0]]:
                scf.contigs.insert(0, contig)
                scf.contigs[0].gap = adj_distance
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


def _make_contigs(perm_container):
    """
    Converts permutations into contigs
    """
    contigs = []
    index = defaultdict(list)
    for perm in perm_container.target_perms_filtered:
        assert len(perm.blocks)

        contigs.append(Contig(perm.chr_name))
        for block in perm.blocks:
            index[block.block_id].append(contigs[-1])
            contigs[-1].blocks.append(block.signed_id())

    return contigs, index


def _calc_n50(scaffolds_lengths, assembly_len):
    n50 = 0
    sum_len = 0
    for l in sorted(scaffolds_lengths, reverse=True):
        sum_len += l
        if sum_len > assembly_len / 2:
            n50 = l
            break
    return n50
