#This module provedes some functions 
#for debug output
#############################################

import os
import shutil

#PUBLIC:
#############################################

#singleton providing global debug configuration
class DebugConfig():
    instance = None

    def __init__(self):
        self.debug_dir = None
        self.debugging = False
        self.gen_to_color = {}
        self.colors = ["red", "green", "blue", "yellow",
                        "cyan", "magnetta"]

    def set_debug_dir(self, debug_dir):
        self.debug_dir = debug_dir
        self.debugging = True
        if os.path.isdir(debug_dir):
            shutil.rmtree(debug_dir)
        os.mkdir(debug_dir)

    def genome_to_color(self, genome_id):
        if genome_id not in self.gen_to_color:
            self.gen_to_color[genome_id] = self.colors[0]
            self.colors = self.colors[1:] + self.colors[:1] #rotate list
        return self.gen_to_color[genome_id]

    @staticmethod
    def get_instance():
        if not DebugConfig.instance:
            DebugConfig.instance = DebugConfig()
        return DebugConfig.instance

    #outputs colored breakpoint graph in "dot" format
    def output_bg_component(self, component, name):
        dot_file = open(os.path.join(self.debug_dir, name), "w")
        dot_file.write("graph {\n")
        for v_1, v_2, edge_data in component.edges(data=True):
            label = ""
            color = "black"
            if "label" in edge_data:
                label = edge_data["label"]
            if "genome_id" in edge_data:
                genome_id = edge_data["genome_id"]
                color = self.genome_to_color(genome_id)

            dot_file.write("{0} -- {1} [color=\"{2}\", label=\"{3}\"];\n"
                                        .format(v_1, v_2, color, label))
        dot_file.write("}\n")
