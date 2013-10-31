from collections import namedtuple, defaultdict
#import sibelia_parser as sp
import networkx as nx
import os
import sys
#import copy


class Permutation:
    def __init__(self, ref_id, chr_id, chr_num, blocks=[]):
        self.ref_id = ref_id
        self.chr_id = chr_id
        self.chr_num = chr_num
        self.blocks = blocks

    def iter_blocks(self, circular=False):
        assert len(self.blocks) > 0
        for block in self.blocks:
            yield block

        if circular:
            yield self.blocks[0]


def find_duplications(ref_perms, target_perms):
    index = defaultdict(set)
    duplications = set()
    for perm in ref_perms + target_perms:
        for block in perm.blocks:
            if perm.ref_id in index[block]:
                duplications.add(block)
            else:
                index[block].add(perm.ref_id)

    return duplications


def filter_perm(perm, to_filter):
    new_perm = Permutation(perm.ref_id, perm.chr_id, perm.chr_num, [])
    for block in perm.blocks:
        if block not in to_filter:
            new_perm.blocks.append(block)
    return new_perm


def parse_blocks_file(ref_id, filename):
    name = ""
    permutations = []
    chr_count = 0
    for line in open(filename, "r").read().splitlines():
        if line.startswith(">"):
            name = line[1:]
        else:
            blocks = line.split(" ")[:-1]
            permutations.append(Permutation(ref_id, name, chr_count, map(int, blocks)))
            chr_count += 1
    return permutations


def parse_config(filename):
    references = {}
    target = {}
    for line in open(filename, "r").read().splitlines():
        if line.startswith("#"):
            continue

        if line.startswith("REF"):
            ref_id, ref_file = line[4:].split("=")
            references[ref_id] = ref_file

        if line.startswith("TARGET"):
            ref_id, ref_file = line[7:].split("=")
            target[ref_id] = ref_file

    return references, target


class PermutationsContainer:
    def __init__(self, config_file):
        self.ref_perms = []
        self.target_perms = []

        ref_files, target_files = parse_config(config_file)
        prefix = os.path.dirname(config_file)
        for ref_id, ref_file in ref_files.iteritems():
            filename = os.path.join(prefix, ref_file)
            self.ref_perms.extend(parse_blocks_file(ref_id, filename))

        for t_id, t_file in target_files.iteritems():
            filename = os.path.join(prefix, t_file)
            self.target_perms.extend(parse_blocks_file(t_id, filename))

        self.duplications = find_duplications(self.ref_perms, self.target_perms)

        self.target_blocks = set()
        for perm in self.target_perms:
            self.target_blocks |= set(perm.blocks)

        to_filter = self.duplications | self.target_blocks
        self.ref_perms_filtered = [filter_perm(p, to_filter) for p in self.ref_perms]
        self.target_perms_filtered = [filter_perm(p, to_filter) for p in self.target_perms]


class BreakpointGraph:
    def __init__(self):
        self.graph = nx.MultiGraph()


    def get_subgraphs(self):
         return nx.connected_component_subgraphs(self.graph)


    def build_from(self, perm_container, circular):
        #TODO: add target here
        for perm in perm_container.ref_perms_filtered:
            prev = None
            for block in perm.iter_blocks(circular):
                if not prev:
                    prev = block
                    continue

                left_block = prev
                right_block = block

                self.graph.add_node(-left_block)
                self.graph.add_node(right_block)
                self.graph.add_edge(-left_block, right_block, ref_id=perm.ref_id)

                prev = block


#TODO: graph output
Colors = ["red", "green", "blue", "yellow", "black"]

def write_colored_dot(graph, dot_file):
    def output_subgraph(subgraph):
        for edge in subgraph.edges(data=True):
            color = Colors[edge[2]["color"] - 1]
            dot_file.write("""{0} -- {1} [color = "{2}"];\n"""
                            .format(edge[0], edge[1], color))

    dot_file.write("graph {\n")
    output_subgraph(graph)
    dot_file.write("}\n")


if __name__ == "__main__":
    pc = PermutationsContainer(sys.argv[1])
    bg = BreakpointGraph()
    bg.build_from(pc, True)
