#This module provides PermutationContainer class
#which stores permutations and provides some filtering
#procedures
######################################################

from collections import defaultdict
import logging
import os

from ragout.shared.debug import DebugConfig

logger = logging.getLogger()
debugger = DebugConfig.get_instance()

class PermException(Exception):
    pass

#PUBLIC:
########################################################

class Permutation:
    def __init__(self, genome_id, chr_id, chr_num, blocks):
        self.genome_id = genome_id
        self.chr_id = chr_id
        self.chr_num = chr_num
        self.blocks = blocks

    def iter_pairs(self):
        for pb, nb in zip(self.blocks[:-1], self.blocks[1:]):
            yield pb, nb


class PermutationContainer:
    #parses permutation files referenced from config and filters duplications
    def __init__(self, permutations_file, config):
        self.ref_perms = []
        self.target_perms = []

        logging.info("Reading permutation file")
        permutations = _parse_blocks_file(permutations_file)
        if not permutations:
            raise PermException("Error reading permutations")

        for p in permutations:
            if p.genome_id not in config.genomes:
                continue
            if p.genome_id not in config.targets:
                self.ref_perms.append(p)
            elif p.genome_id in config.targets:
                self.target_perms.append(p)

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

        self.target_blocks = set()
        for perm in self.target_perms:
            self.target_blocks |= set(map(abs, perm.blocks))

        #filter dupilcated blocks
        self.duplications = _find_duplications(self.ref_perms,
                                               self.target_perms)
        to_hold = self.target_blocks - self.duplications
        self.ref_perms_filtered = [_filter_perm(p, to_hold)
                                      for p in self.ref_perms]
        self.target_perms_filtered = [_filter_perm(p, to_hold)
                                         for p in self.target_perms]
        self.target_perms_filtered = list(filter(lambda p: p.blocks,
                                                 self.target_perms_filtered))

        if debugger.debugging:
            file = os.path.join(debugger.debug_dir, "used_contigs.txt")
            _write_permutations(self.target_perms_filtered, open(file, "w"))


#PRIVATE:
#######################################################

#find duplicated blocks
def _find_duplications(ref_perms, target_perms):
    index = defaultdict(set)
    duplications = set()
    for perm in ref_perms + target_perms:
        for block in map(abs, perm.blocks):
            if perm.genome_id in index[block]:
                duplications.add(block)
            else:
                index[block].add(perm.genome_id)

    return duplications


#filters duplications
def _filter_perm(perm, to_hold):
    new_perm = Permutation(perm.genome_id, perm.chr_id, perm.chr_num, [])
    for block in perm.blocks:
        if abs(block) in to_hold:
            new_perm.blocks.append(block)
    return new_perm


#parses synteny blocks file
def _parse_blocks_file(filename):
    permutations = []
    chr_count = 0
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                tokens = line[1:].split(".", 1)
                if len(tokens) < 2:
                    logger.error("permutation ids in " + filename + " do not "
                                 "follow naming convention: genome.chromosome")
                    return None

                genome_name = tokens[0]
                chr_name = tokens[1]
            else:
                blocks = line.split(" ")[:-1]
                permutations.append(Permutation(genome_name, chr_name,
                                    chr_count, list(map(int, blocks))))
                chr_count += 1
    return permutations


#outputs permutations
def _write_permutations(permutations, out_stream):
    for perm in permutations:
        out_stream.write(">" + perm.chr_id + "\n")
        for block in perm.blocks:
            out_stream.write("{0:+} ".format(block))
        out_stream.write("$\n")
