from collections import namedtuple, defaultdict
import sibelia_parser as sp
import numpy

Colors = ["red", "green", "blue", "yellow", "black"]
Edge = namedtuple("Edge", ["vertex", "color", "distance"])

class Node:
    def __init__(self, n_id, in_assembly):
        self.node_id = n_id
        self.in_assembly = in_assembly
        self.edges = []


def build_graph(sibelia_output):
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

    graph = {}
    for perm in permutations:
        prev = 0
        while abs(perm.blocks[prev]) in duplications:
            prev += 1
        cur = prev + 1
        while True:
            while cur < len(perm.blocks) and abs(perm.blocks[cur]) in duplications:
                cur += 1
            if cur >= len(perm.blocks):
                break

            left_block = perm.blocks[prev]
            right_block = perm.blocks[cur]
            dist = sibelia_output.get_blocks_distance(abs(left_block), abs(right_block), perm.chr_num)

            #if abs(left_block) in contig_index and abs(right_block) in contig_index:
            if -left_block not in graph:
                graph[-left_block] = Node(-left_block, abs(left_block) in contig_index)
            graph[-left_block].edges.append(Edge(right_block, perm.chr_num, dist))

            if right_block not in graph:
                graph[right_block] = Node(right_block, abs(right_block) in contig_index)
            graph[right_block].edges.append(Edge(-left_block, perm.chr_num, dist))

            prev = cur
            cur += 1
    return graph

def output_graph(graph, dot_file, trivial_con=False):
    dot_file.write("graph {\n")
    used_vertexes = set()
    for node in graph.itervalues():
        color = "black" if node.in_assembly else "grey"
        dot_file.write("""{0} [color = "{1}"];\n""".format(node.node_id, color))

    for node_id, node in graph.iteritems():
        for edge in node.edges:
            if edge.vertex not in used_vertexes:
                dot_file.write("""{0} -- {1} [color = "{2}"];\n"""
                                .format(node_id, edge.vertex, Colors[edge.color - 1]))
        used_vertexes.add(node_id)

    if trivial_con:
        for i in xrange(max(graph.keys()) + 1):
            dot_file.write("""{0} -- {1} [color = "black"];\n""".format(i, -i))
    dot_file.write("}")


def get_connected_components(graph):
    con_comp = []

    visited = set()
    def dfs(vertex, component):
        visited.add(vertex)
        component.append(vertex)
        for edge in graph[vertex].edges:
            if edge.vertex not in visited:
                dfs(edge.vertex, component)

    for vertex in graph:
        if vertex not in visited:
            con_comp.append([])
            dfs(vertex, con_comp[-1])

    return con_comp


def median(vals_list):
    return numpy.median(vals_list)


def get_component_of(connected_comps, vertex):
    for con in connected_comps:
        if vertex in con:
            return con
    return None


def vertex_distance(graph, v1, v2):
    edges = filter(lambda e : e.vertex == v2, graph[v1].edges)
    assert edges
    distance = median(map(lambda e: e.distance, edges))
    #print distance
    return int(distance)
