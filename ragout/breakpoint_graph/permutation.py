#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module provides PermutationContainer class
which parses files, stores permutations and
provides some other usefull functions
"""

from collections import defaultdict
import logging
import os
import math
from copy import deepcopy
from itertools import chain

from ragout.shared.debug import DebugConfig
from ragout.shared import config
from ragout.shared.datatypes import Block, Permutation
import ragout.breakpoint_graph.repeat_resolver as rr

logger = logging.getLogger()
debugger = DebugConfig.get_instance()


class PermException(Exception):
    pass


class PermutationContainer:
    def __init__(self, block_coords_file, recipe,
                 resolve_repeats, conservative, phylogeny):
        """
        Parses permutation files referenced from recipe and filters duplications
        """
        self.ref_perms = []
        self.target_perms = []
        self.recipe = recipe

        logging.debug("Reading permutation file")
        permutations = _parse_blocks_coords(block_coords_file)
        if not permutations:
            raise PermException("Error reading permutations")

        has_sequences = set()
        draft_names = set()
        for p in permutations:
            if p.genome_name not in recipe["genomes"]:
                continue
            p.draft = recipe["genomes"][p.genome_name]["draft"]
            if p.draft:
                draft_names.add(p.genome_name)
            p.circular = recipe["genomes"][p.genome_name]["circular"]

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

        logger.debug("Read {0} reference sequences"
                     .format(len(self.ref_perms)))
        if not len(self.ref_perms):
            raise PermException("No synteny blocks found in "
                                "reference sequences")
        logger.debug("Read {0} target sequences"
                     .format(len(self.target_perms)))
        if not len(self.target_perms):
            raise PermException("No synteny blocks found in "
                                "target sequences")

        self.filter_indels(not conservative)
        logger.debug("{0} target sequences left after indel filtering"
                                        .format(len(self.target_perms)))

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
        logger.debug("{0} target sequences left after repeat filtering"
                     .format(len(self.target_perms)))

        if debugger.debugging:
            file = os.path.join(debugger.debug_dir, "used_contigs.txt")
            _write_permutations(self.target_perms, open(file, "w"))


    def filter_indels(self, allow_ref_indels):
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
            to_keep = filter(lambda b: multiplicity[b] >= num_genomes,
                             multiplicity)

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


def _parse_permutations(filename):
    """
    Parses a file with signed permutations
    """
    permutations = []
    chr_count = 0
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                tokens = line[1:].split(".", 1)
                if len(tokens) != 2:
                    raise PermException("permutation ids in " + filename +
                                        " do not follow naming convention: " +
                                        "'genome.chromosome'")

                genome_name, chr_name = tokens
            else:
                permutations.append(Permutation(genome_name, chr_name,
                                    chr_count, None, []))
                blocks_ids = map(int, line.split(" ")[:-1])
                for b in blocks_ids:
                    permutations[-1].blocks.append(Block(abs(b),
                                                         math.copysign(1, b)))
                chr_count += 1
    return permutations


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
                                                 chr_id, int(chr_size), [])

            else:
                if line.startswith("Seq_id") or line.startswith("-"):
                    continue

                if line.startswith("Block"):
                    block_id = int(line.split(" ")[1][1:])
                    continue

                seq_id, sign, start, end, length = line.split("\t")
                if sign == "-":
                    start, end = end, start
                if int(end) < int(start):
                    raise PermException("Error in permutations file format")

                sign_num = 1 if sign == "+" else -1
                perm_by_id[seq_id].blocks.append(Block(block_id, sign_num,
                                                      int(start), int(end)))

    for perm in perm_by_id.values():
        perm.blocks.sort(key=lambda b: b.start)

    out_perms = list(filter(lambda b: len(b.blocks), perm_by_id.values()))
    return out_perms


def _check_coverage(permutations):
    """
    Checks if synteny blocks coverage is acceptable
    """
    by_genome = defaultdict(list)
    for perm in permutations:
        by_genome[perm.genome_name].append(perm)

    for genome_name, genome_perms in by_genome.items():
        total_length = 0
        total_covered = 0
        for perm in genome_perms:
            total_length += perm.chr_len
            for block in perm.blocks:
                total_covered += block.length()

        coverage = float(total_covered) / total_length
        logger.debug("\"{0}\" synteny blocks coverage: {1:2.4}%"
                     .format(genome_name, 100 * coverage))
        if coverage < config.vals["min_synteny_coverage"]:
            logger.warning("\"{0}\" synteny blocks coverage ({1:2.4}%) "
                           "is too low -- results can be inaccurate. "
                           "Possible reasons are: \n\n"
                           "1. Genome is too distant from others\n"
                           "2. Synteny block parameters are chosen incorrectly"
                           "\n\nTry to change synteny blocks paramsters or "
                           "remove this genome from the comparison"
                           .format(genome_name, 100 * coverage))


def _write_permutations(permutations, out_stream):
    """
    Outputs permutations
    """
    for perm in permutations:
        out_stream.write(">" + perm.chr_name + "\n")
        for block in perm.blocks:
            out_stream.write("{0:+} ".format(block.signed_id()))
        out_stream.write("$\n")
