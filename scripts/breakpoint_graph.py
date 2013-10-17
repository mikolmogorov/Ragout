from collections import namedtuple, defaultdict
import sibelia_parser as sp
import networkx as nx
import copy

MIN_DEGREE = 1

class BreakpointGraph:
    def __init__(self):
        self.graph = nx.MultiGraph()


    def get_subgraphs(self):
         return nx.connected_component_subgraphs(self.graph)


    def get_subgraph_with(self, node_id):
        assert self.unresolved_subgraphs
        for subgr in self.unresolved_subgraphs:
            if subgr.has_node(node_id):
                return subgr
        return None


    def build_from(self, sibelia_output):
        blocks_coords = sibelia_output.blocks_info
        contig_index = sp.build_contig_index(sibelia_output.contigs)

        duplications = sibelia_output.get_duplications()

        circular_perm = copy.deepcopy(sibelia_output.permutations)
        for perm in circular_perm:
            perm = copy.copy(perm)
            perm.blocks.append(perm.blocks[0])

        for perm in circular_perm:
            prev = 0
            while (abs(perm.blocks[prev]) in duplications
                    or not abs(perm.blocks[prev]) in contig_index):
                #print prev
                prev += 1
            cur = prev + 1
            while True:

                end = False
                while True:
                    if cur >= len(perm.blocks):
                        end = True
                        break
                    in_assembly = abs(perm.blocks[cur]) in contig_index
                    duplication = abs(perm.blocks[cur]) in duplications
                    #degree = block_count[abs(perm.blocks[cur])]
                    #if not duplication and in_assembly and degree >= MIN_DEGREE:
                    if not duplication and in_assembly:
                        break
                    cur += 1

                if end:
                    break

                left_block = perm.blocks[prev]
                right_block = perm.blocks[cur]
                dist = sibelia_output.get_blocks_distance(abs(left_block), abs(right_block), perm.chr_num)
                ref_id = sibelia_output.chr_id_to_ref_id[perm.chr_id]

                self.graph.add_node(-left_block)
                self.graph.add_node(right_block)
                self.graph.add_edge(-left_block, right_block, color=perm.chr_num,
                                    distance=dist, ref_id=ref_id)

                prev = cur
                cur += 1


    def vertex_distance(self, v1, v2):
        edges = self.graph[v1][v2].values()
        distance = weighted_distance(edges)
        return int(distance)



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


def weighted_distance(edges):
    assert edges
    return sorted(edges, key=lambda e: e["color"])[0]["distance"]
