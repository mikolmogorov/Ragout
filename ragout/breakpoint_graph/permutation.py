#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module provides PermutationContainer class
which parses files, stores permutations and
provides some other usefull functions
"""

from __future__ import absolute_import
from __future__ import division
from collections import defaultdict
import logging
import os
from copy import deepcopy

from ragout.shared.debug import DebugConfig
from ragout.shared import config
from ragout.shared.datatypes import Block, Permutation, output_permutations
import ragout.breakpoint_graph.repeat_resolver as rr
from ragout.six.moves import filter

logger = logging.getLogger()
debugger = DebugConfig.get_instance()


class PermException(Exception):
    pass


class PermutationContainer:
    def __init__(self, block_coords_file, recipe,
                 resolve_repeats, allow_ref_indels, phylogeny):
        """
        Parses permutation files referenced from recipe and filters duplications
        """
        self.ref_perms = []
        self.target_perms = []
        self.recipe = recipe

        logging.info("Reading %s", block_coords_file)
        permutations = _parse_blocks_coords(block_coords_file)

        has_sequences = set()
        draft_names = set()
        for p in permutations:
            if p.genome_name not in recipe["genomes"]:
                continue
            p.draft = recipe["genomes"][p.genome_name]["draft"]
            if p.draft:
                draft_names.add(p.genome_name)

            has_sequences.add(p.genome_name)
            if p.genome_name == recipe["target"]:
                self.target_perms.append(p)
            else:
                self.ref_perms.append(p)

        for genome in recipe["genomes"]:
            if genome not in has_sequences:
                raise PermException("No sequences read for genome {0}. Check "
                                    "recipe for correctness.".format(genome))
        _check_coverage(self.ref_perms + self.target_perms)

        logger.debug("Read %d reference sequences", len(self.ref_perms))
        if not len(self.ref_perms):
            raise PermException("No synteny blocks found in "
                                "reference sequences")
        logger.debug("Read %d target sequences", len(self.target_perms))
        if not len(self.target_perms):
            raise PermException("No synteny blocks found in "
                                "target sequences")

        self._filter_indels(allow_ref_indels)
        logger.debug("%d target sequences left after indel filtering",
                     len(self.target_perms))

        repeats = _find_repeats(self.ref_perms + self.target_perms)
        ###
        if resolve_repeats:
            if phylogeny is None:
                raise PermException("Resolving repeats with "
                                    "yet unknown phylogeny")
            rr.resolve_repeats(self.ref_perms, self.target_perms,
                               repeats, phylogeny, draft_names)
        ###
        self._filter_repeats(repeats)
        logger.debug("%d target sequences left after repeat filtering",
                     len(self.target_perms))
        if not len(self.target_perms):
            raise PermException("No synteny blocks found in the target "
                                "genome after repeat/indel filtering.")

        if debugger.debugging:
            filtered_file = os.path.join(debugger.debug_dir, "filtered_contigs.txt")
            output_permutations(self.target_perms, filtered_file)

    def _filter_indels(self, allow_ref_indels):
        """
        Keep only blocks that appear in target and
        all references (or one reference, if allow_ref_indels is set)
        """
        multiplicity = defaultdict(int)
        target_blocks = set()
        reference_blocks = set()
        def process(perms, block_set):
            for perm in perms:
                for block in perm.blocks:
                    multiplicity[block.block_id] += 1
                    block_set.add(block.block_id)

        process(self.target_perms, target_blocks)
        process(self.ref_perms, reference_blocks)

        if allow_ref_indels:
            to_keep = target_blocks.intersection(reference_blocks)
        else:
            num_genomes = len(self.recipe["genomes"])
            #to_keep = set(filter(lambda b: multiplicity[b] >= num_genomes,
            #                     multiplicity))
            to_keep = set([b for b in multiplicity if multiplicity[b] >= num_genomes])

        self.ref_perms = _filter_permutations(self.ref_perms, to_keep)
        self.target_perms = _filter_permutations(self.target_perms, to_keep)


    def _filter_repeats(self, repeats):
        """
        Filters repetitive blocks
        """
        self.target_perms = _filter_permutations(self.target_perms, repeats,
                                                 inverse=True)
        self.ref_perms = _filter_permutations(self.ref_perms, repeats,
                                              inverse=True)


def _find_repeats(permutations):
    """
    Returns a set of repetitive blocks
    """
    index = defaultdict(set)
    repeats = set()
    for perm in permutations:
        for block in perm.blocks:
            if perm.genome_name in index[block.block_id]:
                repeats.add(block.block_id)
            else:
                index[block.block_id].add(perm.genome_name)
    return repeats


def _filter_permutations(permutations, blocks, inverse=False):
    """
    Filters given blocks from permutations
    """
    new_perms = []
    filter_func = lambda b: (b.block_id in blocks) != inverse

    for perm in permutations:
        new_blocks = list(filter(filter_func, perm.blocks))
        if new_blocks:
            new_perms.append(deepcopy(perm))
            new_perms[-1].blocks = new_blocks
    return new_perms


def _parse_blocks_coords(filename):
    """
    Parses a file with blocks coords
    """
    perm_by_id = {}
    with open(filename, "r") as f:
        header = True
        for line in f:
            line = line.strip()
            if not line:
                continue

            if header:
                if line.startswith("Seq_id"):
                    continue

                if line.startswith("-"):
                    header = False
                    continue

                chr_id, chr_size, seq_name = line.split("\t")
                tokens = seq_name.split(".", 1)
                if len(tokens) != 2:
                    raise PermException("permutation ids in " + filename +
                                        " do not follow naming convention: " +
                                        "'genome.chromosome'")

                genome_name, chr_name = tokens
                perm_by_id[chr_id] = Permutation(genome_name, chr_name,
                                                 int(chr_size), [])

            else:
                if line.startswith("Seq_id") or line.startswith("-"):
                    continue

                if line.startswith("Block"):
                    block_id = int(line.split(" ")[1][1:])
                    continue

                seq_id, sign, start, end, _length = line.split("\t")
                if sign == "-":
                    start, end = end, start
                if int(end) < int(start):
                    raise PermException("Error in permutations file format")

                sign_num = 1 if sign == "+" else -1
                perm_by_id[seq_id].blocks.append(Block(block_id, sign_num,
                                                      int(start), int(end)))

    for perm in perm_by_id.values():
        perm.blocks.sort(key=lambda b: b.start)

    #out_perms = list(filter(lambda b: len(b.blocks), perm_by_id.values()))
    out_perms = [b for b in  perm_by_id.values() if len(b.blocks)]
    if not len(out_perms):
        raise PermException("Permutations file is empty")
    return out_perms


def _check_coverage(permutations):
    """
    Checks if synteny blocks coverage is acceptable
    """
    by_genome = defaultdict(list)
    for perm in permutations:
        by_genome[perm.genome_name].append(perm)

    warning_genomes = []
    for genome_name, genome_perms in by_genome.items():
        total_length = 0
        total_covered = 0
        for perm in genome_perms:
            total_length += perm.length()
            for block in perm.blocks:
                total_covered += block.length()

        coverage = float(total_covered) / total_length
        logger.info("\"{0}\" synteny blocks coverage: {1:2.4}%"
                     .format(genome_name, 100 * coverage))
        if coverage < config.vals["min_synteny_coverage"]:
            warning_genomes.append((genome_name, coverage))

    if warning_genomes:
        #for genome_name, coverage in warning_genomes:
        #    logger.warning("\"{0}\" synteny blocks coverage: {1:2.4}%"
        #                 .format(genome_name, 100 * coverage))
        logger.warning("Some genomes have low synteny block coverage "
                       "- results can be inaccurate. "
                       "Possible reasons: \n"
                       "\t1. Genome is too distant from others\n"
                       "\t2. Synteny block parameters are chosen incorrectly"
                       "\n\tTry to change synteny blocks paramsters or "
                       "remove this genome from the comparison")
