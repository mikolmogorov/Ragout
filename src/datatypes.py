#This module provides some common data structures
#################################################

from collections import namedtuple

#PUBLIC:
#################################################

class Scaffold:
    def __init__(self, name):
        self.left = 0
        self.right = 0
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
    def __init__(self, name, sign=1, gap=0):
        self.name = name
        self.sign = sign
        self.gap = gap
        self.blocks = []

    @staticmethod
    def from_sting(string):
        return Contig(*parse_contig_name(string))

    def __str__(self):
        sign = "+" if self.sign > 0 else "-"
        return sign + self.name

def parse_contig_name(string):
    if string[0] not in ["+", "-"]:
        return None

    sign = 1 if string[0] == "+" else -1
    return string[1:], sign
