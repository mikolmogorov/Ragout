from collections import namedtuple, defaultdict
import sibelia_parser as sp
import networkx as nx

Colors = ["red", "green", "blue", "yellow", "black"]


class BreakpointGraph:
    def __init__(self):
        self.graph = nx.MultiGraph()
        self.unresolved_subgraphs = None


    def mark_resolved(self, node_ids):
        new = filter(lambda x: x[0] not in node_ids, enumerate(self.unresolved_subgraphs))
        self.unresolved_subgraphs = map(lambda x: x[1], new)


    def get_subgraph_with(self, node_id):
        assert self.unresolved_subgraphs
        for subgr in self.unresolved_subgraphs:
            if subgr.has_node(node_id):
                return subgr
        return None


    def build_from(self, sibelia_output):
        permutations = sibelia_output.permutations
        blocks_coords = sibelia_output.blocks_info
        contig_index = sibelia_output.build_contig_index()

        #find duplications
        duplications = set()
        for perm in permutations:
            current = set()
            for block in perm.blocks:
                if abs(block) in current:
                    duplications.add(abs(block))
                current.add(abs(block))
        print "Duplications found: ", duplications

        for perm in permutations:
            prev = 0
            while abs(perm.blocks[prev]) in duplications:
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
                    #if not duplication and in_assembly:
                    if not duplication:
                        break
                    cur += 1

                if end:
                    break

                left_block = perm.blocks[prev]
                right_block = perm.blocks[cur]
                dist = sibelia_output.get_blocks_distance(abs(left_block), abs(right_block), perm.chr_num)

                self.graph.add_node(-left_block, in_assembly=(abs(left_block) in contig_index))
                self.graph.add_node(right_block, in_assembly=(abs(right_block) in contig_index))
                self.graph.add_edge(-left_block, right_block, color=perm.chr_num, distance=dist)

                prev = cur
                cur += 1


    def write_dot(self, dot_file):
        def output_subgraph(subgraph):
            for edge in subgraph.edges(data=True):
                color = Colors[edge[2]["color"] - 1]
                dot_file.write("""{0} -- {1} [color = "{2}"];\n"""
                                .format(edge[0], edge[1], color))

        dot_file.write("graph {\n")
        for node in self.graph.nodes(data=True):
            color = "black" if node[1]["in_assembly"] else "grey"
            dot_file.write("""{0} [color = "{1}"];\n""".format(node[0], color))

        unresolved_nodes = map(lambda s: s.nodes(), self.unresolved_subgraphs)
        unresolved_nodes = sum(unresolved_nodes, [])
        resolved_nodes = filter(lambda n: n not in unresolved_nodes, self.graph.nodes())
        resolved_subgr = self.graph.subgraph(resolved_nodes)

        dot_file.write("""subgraph cluster_0 {\nlabel = "resolved components"\n""")
        output_subgraph(resolved_subgr)

        dot_file.write("""}\nsubgraph cluster_1 {\nlabel = "unresolved components"\n""")
        for subgr in self.unresolved_subgraphs:
            output_subgraph(subgr)
        dot_file.write("}\n}\n")


    def in_assembly(self, node_id):
        return self.graph.node[node_id]["in_assembly"]


    def vertex_distance(self, v1, v2):
        edges = self.graph[v1][v2].values()
        distance = weighted_distance(edges)
        return int(distance)


    def get_unresolved_subgraphs(self):
        if not self.unresolved_subgraphs:
            self.unresolved_subgraphs = nx.connected_component_subgraphs(self.graph)
        return self.unresolved_subgraphs


def weighted_distance(edges):
    assert edges
    return sorted(edges, key=lambda e: e["color"])[0]["distance"]
