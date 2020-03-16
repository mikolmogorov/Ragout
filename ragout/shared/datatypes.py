#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module provides some common data structures
"""

from __future__ import absolute_import
from __future__ import division
from copy import copy
from ragout.six.moves import zip


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
    def __init__(self, genome_name, chr_name, seq_len, blocks):
        self.genome_name = genome_name
        self.chr_name = chr_name
        self.seq_start = 0
        self.seq_end = seq_len
        self.seq_len = seq_len
        self.blocks = blocks
        self.repeat_id = 0
        self.draft = False

    def length(self):
        assert self.seq_end > self.seq_start
        return self.seq_end - self.seq_start

    def name(self):
        if self.seq_start == 0 and self.seq_end == self.seq_len:
            return self.chr_name
        else:
            return "{0}[{1}:{2}]".format(self.chr_name, self.seq_start,
                                         self.seq_end)

    def iter_pairs(self):
        for pb, nb in zip(self.blocks[:-1], self.blocks[1:]):
            yield pb, nb

    def __lt__(self, other):
        return repr(self) < repr(other)

    def __repr__(self):
        return ("[{0}, {1}, {2}, b:{3}, e:{4}]"
                    .format(self.genome_name, self.chr_name,
                            [b.signed_id() for b in self.blocks],
                            self.seq_start, self.seq_end))


def output_permutations(permutations, out_file):
    with open(out_file, "w") as f:
        for perm in permutations:
            f.write(">" + perm.name() + "\n")
            for block in perm.blocks:
                f.write("{0:+} ".format(block.signed_id()))
            f.write("$\n")


class Contig:
    def __init__(self, permutation, sign, link, _dummy_param):
        self.perm = permutation
        self.sign = sign
        if link is None:
            link = Link(0, [])
        self.link = link

    @staticmethod
    def with_perm(permutation, sign=1, link=None):
        return Contig(permutation, sign, link, None)

    @staticmethod
    def with_sequence(seq_name, seq_len, sign=1, link=None):
        dummy_perm = Permutation(None, seq_name, seq_len, None)
        return Contig.with_perm(dummy_perm, sign, link)

    def left_end(self):
        return (self.perm.blocks[0].signed_id() if self.sign > 0
                else -self.perm.blocks[-1].signed_id())

    def right_end(self):
        return (-self.perm.blocks[-1].signed_id() if self.sign > 0
                else self.perm.blocks[0].signed_id())

    def left_gap(self):
        return (self.perm.blocks[0].start if self.sign > 0
                else self.perm.length() - self.perm.blocks[-1].end)

    def right_gap(self):
        return (self.perm.length() - self.perm.blocks[-1].end
                if self.sign > 0 else self.perm.blocks[0].start)

    def reverse_copy(self):
        contig = copy(self)
        contig.sign = -contig.sign
        return contig

    def signed_perm(self):
        if self.sign > 0:
            return [b.signed_id() for b in self.perm.blocks]
        else:
            return [-b.signed_id() for b in self.perm.blocks[::-1]]

    def name(self):
        return self.perm.name()

    def signed_name(self):
        sign = "+" if self.sign > 0 else "-"
        return sign + self.name()

    def name_with_coords(self):
        return self.perm.chr_name, self.perm.seq_start, self.perm.seq_end

    def trim_left(self, trim_len):
        if self.sign > 0:
            self.perm.seq_start += trim_len
        else:
            self.perm.seq_end -= trim_len

    def trim_right(self, trim_len):
        if self.sign > 0:
            self.perm.seq_end -= trim_len
        else:
            self.perm.seq_start += trim_len

    def length(self):
        return self.perm.length()


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


def output_scaffolds_premutations(scaffolds, out_file):
    with open(out_file, "w") as f:
        for scf in scaffolds:
            blocks = []
            for contig in scf.contigs:
                blocks.extend(contig.signed_perm())

            f.write(">" + scf.name + "\n")
            for block in blocks:
                f.write("{0:+} ".format(block))
            f.write("$\n")
