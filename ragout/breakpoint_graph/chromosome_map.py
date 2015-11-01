#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module uses chromosome map to distinguish
between true and chimeric adjacencies
"""

from __future__ import print_function
from collections import namedtuple, defaultdict
import logging
from copy import copy, deepcopy

from ragout.breakpoint_graph.chimera_detector import _break_permutation

ChrBreak = namedtuple("ChrBreak", ["chr_id", "pos"])
ChrJoin = namedtuple("ChrJoin", ["chr_id_1", "pos_1", "chr_id_2", "pos_2"])

logger = logging.getLogger()

class ChromosomeMap(object):
    def __init__(self, map_file):
        self.breaks = []
        self.joins = []
        self.ref_name = None
        self._read_map(map_file)

    def _read_map(self, map_file):
        for line in open(map_file, "r"):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith(".reference"):
                self.ref_name = line.split("=")[1].strip()

            elif line.startswith(".break"):
                val = line.split("=")[1].strip()
                chr_id, position = val.split(":")
                self.breaks.append(ChrBreak(chr_id, _pos_num(position)))

            elif line.startswith(".join"):
                val = line.split("=")[1].strip()
                left, right = val.split("-")
                chr_1, pos_1 = left.split(":")
                chr_2, pos_2 = right.split(":")
                self.joins.append(ChrJoin(chr_1, _pos_num(pos_1),
                                          chr_2, _pos_num(pos_2)))

    def fix_container_and_graph(self, perm_container, bp_graph):
        breaks_by_chr = defaultdict(list)
        for br in self.breaks:
            breaks_by_chr[br.chr_id].append(br.pos)

        breakpoints = []
        for perm in perm_container.ref_perms:
            if perm.genome_name == self.ref_name:
                for br in breaks_by_chr[perm.chr_name]:
                    breakpoints.append(_get_breakpoint(perm, br))

        target_breaks = defaultdict(list)
        for node_1, node_2 in breakpoints:
            for edge in bp_graph.bp_graph[node_1][node_2].values():
                if edge["genome_id"] == bp_graph.target:
                    target_breaks[edge["chr_name"]].append(edge["start"])
            bp_graph.bp_graph.remove_node(node_1)
            bp_graph.bp_graph.remove_node(node_2)

        new_target_perms = []
        for perm in perm_container.target_perms:
            if perm.chr_name in target_breaks:
                new_target_perms.extend(_break_permutation(perm,
                                                target_breaks[perm.chr_name]))
                logger.debug("Extra break")
            else:
                new_target_perms.append(perm)
        perm_container.target_perms = new_target_perms


def _get_breakpoint(perm, position):
    for prev_block, next_block in zip(perm.blocks[:-1], perm.blocks[1:]):
        if prev_block.start <= position < next_block.end:
            return (-prev_block.signed_id(), next_block.signed_id())
    return None


def _pos_num(str_pos):
    if str_pos == "start":
        return 0
    elif str_pos == "end":
        return -1
    else:
        return int(str_pos)
