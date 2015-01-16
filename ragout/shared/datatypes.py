#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module provides some common data structures
"""

from collections import namedtuple
from copy import deepcopy


class Block:
    """
    Represents synteny block
    """
    def __init__(self, block_id, sign, start=None, end=None):
        self.block_id = block_id
        self.sign = sign
        self.start = start
        self.end = end

    def length(self):
        if self.start is None or self.end is None:
            return None

        assert self.end >= self.start
        return self.end - self.start

    def signed_id(self):
        return self.block_id * self.sign


class Permutation:
    """
    Represents signed permutation
    """
    def __init__(self, genome_name, chr_name, chr_id, chr_len, blocks):
        self.genome_name = genome_name
        self.chr_name = chr_name
        self.chr_id = chr_id
        self.chr_len = chr_len
        self.blocks = blocks
        #self.chr_index = None

    def iter_pairs(self):
        for pb, nb in zip(self.blocks[:-1], self.blocks[1:]):
            yield pb, nb


class Contig:
    """
    Contig data structure for more convenient use
    """
    def __init__(self, seq_name, permutation, sign=1):
        self.perm = permutation
        self.seq_name = seq_name
        self.sign = sign
        self.link = Link(0, [])

    def left(self):
        return self.perm.blocks[0].signed_id()

    def right(self):
        return self.perm.blocks[-1].signed_id()

    def reverse_copy(self):
        contig = deepcopy(self)
        contig.sign = -contig.sign
        return contig

    def signed_perm(self):
        if self.sign > 0:
            return list(map(lambda b: b.signed_id(), self.perm.blocks))
        else:
            return list(map(lambda b: -b.signed_id(), self.perm.blocks[::-1]))

    def __str__(self):
        sign = "+" if self.sign > 0 else "-"
        return sign + self.seq_name


class Link:
    """
    Represens an adjancency between teo contigs
    """
    def __init__(self, gap, supporting_genomes):
        self.gap = gap
        self.supporting_genomes = supporting_genomes
        self.supporting_assembly = False


class Scaffold:
    def __init__(self, name):
        self.left = None
        self.right = None
        self.contigs = []
        self.name = name

    @staticmethod
    def with_contigs(name, left, right, contigs):
        scf = Scaffold(name)
        scf.left = left
        scf.right = right
        scf.contigs = contigs
        return scf
