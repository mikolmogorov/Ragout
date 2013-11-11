import os
import shutil

class DebugConfig():
    instance = None

    def __init__(self):
        self.debug_dir = None
        self.debugging = False

    def set_debug_dir(self, debug_dir):
        self.debug_dir = debug_dir
        self.debugging = True
        if os.path.isdir(debug_dir):
            shutil.rmtree(debug_dir)
        os.mkdir(debug_dir)

    @staticmethod
    def get_writer():
        if not DebugConfig.instance:
            DebugConfig.instance = DebugConfig()
        return DebugConfig.instance


def write_dot(graph, dot_file):
    colors = ["red", "green", "blue", "yellow", "black"]
    ref_to_color = {}

    dot_file.write("graph {\n")
    for v_1, v_2, edge_data in graph.edges(data=True):
        label = ""
        color = "black"
        if "label" in edge_data:
            label = edge_data["label"]
        if "ref_id" in edge_data:
            ref_id = edge_data["ref_id"]
            if ref_id not in ref_to_color:
                ref_to_color[ref_id] = colors[0]
                del colors[0]
            color = ref_to_color[ref_id]

        dot_file.write("{0} -- {1} [color=\"{2}\", label=\"{3}\"];\n"
                            .format(v_1, v_2, color, label))
    dot_file.write("}\n")
