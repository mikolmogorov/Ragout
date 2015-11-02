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
ChrJoin = namedtuple("ChrJoin", ["chr_left", "chr_right"])

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
                self.joins.append(ChrJoin(left, right))

    def fix_container_and_graph(self, perm_container, bp_graph):
        breaks_by_chr = defaultdict(list)
        for br in self.breaks:
            breaks_by_chr[br.chr_id].append(br.pos)

        breakpoints = []
        endpoints = {}
        for perm in perm_container.ref_perms:
            if perm.genome_name == self.ref_name:
                bp_list = [perm.blocks[0].signed_id()]
                for br in sorted(breaks_by_chr[perm.chr_name]):
                    breakpoints.append(_get_breakpoint(perm, br, bp_graph))
                    bp_list.append(breakpoints[-1][0])
                    bp_list.append(breakpoints[-1][1])
                bp_list.append(-perm.blocks[-1].signed_id())

                if len(bp_list) == 2:
                    endpoints[perm.chr_name] = tuple(bp_list)
                else:
                    for frag_num in xrange(len(bp_list) / 2):
                        frag_name = "{0}.{1}".format(perm.chr_name, frag_num + 1)
                        endpoints[frag_name] = tuple(bp_list[frag_num * 2 :
                                                             frag_num * 2 + 2])

        #apply breaks to breakpoint graph
        target_breaks = defaultdict(list)
        for node_1, node_2 in breakpoints:
            for edge in bp_graph.bp_graph[node_1][node_2].values():
                if edge["genome_id"] == bp_graph.target:
                    target_breaks[edge["chr_name"]].append(edge["start"])
            _disconnect(bp_graph.bp_graph, node_1)
            _disconnect(bp_graph.bp_graph, node_2)
            bp_graph.bp_graph.add_edge(node_1, node_2, chr_name=None,
                                       genome_id=self.ref_name, infinity=True)

        #apply breaks to target permutations
        #new_target_perms = []
        #for perm in perm_container.target_perms:
        #    if perm.chr_name in target_breaks:
        #        new_target_perms.extend(_break_permutation(perm,
        #                                        target_breaks[perm.chr_name]))
        #        logger.debug("Extra break")
        #    else:
        #        new_target_perms.append(perm)
        #perm_container.target_perms = new_target_perms

        #apply joins to breakpoint graph
        for join in self.joins:
            node_1 = endpoints[join.chr_left][1]
            node_2 = endpoints[join.chr_right][0]
            _disconnect(bp_graph.bp_graph, node_1)
            _disconnect(bp_graph.bp_graph, node_2)
            bp_graph.bp_graph.add_edge(node_1, node_2, chr_name="chr_map",
                                       start=0, end=1000000,
                                       genome_id=self.ref_name, infinity=False)


def _disconnect(graph, node):
   for neighbor in graph.neighbors(node):
       graph.remove_edges_from([(node, neighbor)] * 10)


def _get_breakpoint(perm, position, bp_graph):
    breaking = False
    for prev_block, next_block in zip(perm.blocks[:-1], perm.blocks[1:]):
        if prev_block.start <= position < next_block.end:
            breaking = True

        if breaking:
            v1, v2 = -prev_block.signed_id(), next_block.signed_id()
            genome_ids = set(bp_graph.genomes_support(v1, v2))
            if bp_graph.target not in genome_ids:
                return (v1, v2)
    return None


def _pos_num(str_pos):
    if str_pos == "start":
        return 0
    elif str_pos == "end":
        return -1
    else:
        return int(str_pos)
