#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module tries to detect chimerias in
input scaffolds
"""

from __future__ import print_function

from ragout.breakpoint_graph.breakpoint_graph import BreakpointGraph

class ChimeraDetector(object):
    def __init__(self, breakpoint_graph):
        self.graph = breakpoint_graph

    def get_chimeric_adj(self):
        chimeric_adj = []

        subgraphs = self.graph.connected_components()
        for subgr in subgraphs:
            #if len(subgr.bp_graph) > 100:
            #    continue
            for (u, v) in subgr.bp_graph.edges_iter():
                genomes = subgr.supporting_genomes(u, v)
                if len(genomes) == 1 and genomes[0] == self.graph.target:
                    print(len(subgr.bp_graph))
                    if not subgr.alternating_cycle(u, v):
                        chimeric_adj.append((u, v))
                        print("chimeric", u, v)
                    else:
                        print("non chimeric", u, v)

        return chimeric_adj
