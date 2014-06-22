#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module provides some common data structures
"""

from collections import namedtuple
from copy import copy


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

class Contig:
    def __init__(self, name, sign=1):
        self.name = name
        self.sign = sign
        self.gap = 0
        self.blocks = []

    def reverse(self):
        contig = copy(self)
        contig.sign = -contig.sign
        return contig

    def __str__(self):
        sign = "+" if self.sign > 0 else "-"
        return sign + self.name
