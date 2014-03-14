#This module provedes some functions 
#for debug output
#############################################

import os
import shutil
import networkx as nx
from Bio import Phylo

#PUBLIC:
#############################################

#singleton providing global debug configuration
class DebugConfig():
    instance = None

    def __init__(self):
        self.debug_dir = None
        self.debugging = False
        #self.gen_to_color = {}
        #self.colors = ["red", "green", "blue", "yellow",
        #               "cyan", "magnetta"]

    #also enables debugging
    def set_debug_dir(self, debug_dir):
        #try:
        #    import pygraphviz
        #    import pylab
        #except ImportError:
        #    raise Exception("Debugging requires pygraphviz and matplotlib")

        self.debug_dir = debug_dir
        self.debugging = True
        if os.path.isdir(debug_dir):
            shutil.rmtree(debug_dir)
        os.mkdir(debug_dir)

    """
    def genome_to_color(self, genome_id):
        if genome_id not in self.gen_to_color:
            self.gen_to_color[genome_id] = self.colors[0]
            self.colors = self.colors[1:] + self.colors[:1] #rotate list
        return self.gen_to_color[genome_id]
    """

    @staticmethod
    def get_instance():
        if not DebugConfig.instance:
            DebugConfig.instance = DebugConfig()
        return DebugConfig.instance

    #outputs colored breakpoint graph
    """
    def draw_breakpoint_graph(self, graph, name, weights={}):
        graph_to_draw = nx.MultiGraph()
        for v_1, v_2, edge_data in graph.edges_iter(data=True):
            assert "genome_id" in edge_data
            color = self.genome_to_color(edge_data["genome_id"])
            graph_to_draw.add_nodes_from([v_1, v_2])
            graph_to_draw.add_edge(v_1, v_2, color=color)
        for (v1, v2), weight in weights.iteritems():
            label = "{0:4.2f}".format(weight)
            graph_to_draw[v1][v2][0]["label"] = label

        out_file = os.path.join(self.debug_dir, name)
        agraph = nx.to_agraph(graph_to_draw)
        agraph.layout()
        agraph.draw(out_file)
    """

    """
    def output_phylogeny(self, phylogeny, name):
        import pylab
        for clade in phylogeny.tree.find_clades():
            if clade.is_terminal():
                clade.color = self.genome_to_color(clade.name)
        phylogeny.tree.ladderize()
        pylab.rcParams["lines.linewidth"] = 3.0
        Phylo.draw(phylogeny.tree, do_show=False)

        out_file = os.path.join(self.debug_dir, name)
        pylab.savefig(out_file)
    """
