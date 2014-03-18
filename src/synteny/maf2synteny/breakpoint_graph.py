from __future__ import print_function
from collections import defaultdict
from itertools import izip
import sys

from permutations import Block


class Edge:
    def __init__(self, left_node, right_node, seq_id):
        self.seq_id = seq_id
        self.left_node = left_node
        self.right_node = right_node
        self.left_pos = 0
        self.right_pos = 0
        self.prev_edge = None
        self.next_edge = None

    def has_node(self, node):
        return self.left_node == node or self.right_node == node

    def __str__(self):
        left = str(self.left_node) if self.left_node != sys.maxint else "inf"
        right = str(self.right_node) if self.right_node != sys.maxint else "inf"
        return "({0}, {1}, {2})".format(left, right, self.seq_id)


class Node:
    def __init__(self):
        self.edges = []


class BreakpointGraph:
    def __init__(self):
        self.origins = {}
        self.nodes = defaultdict(Node)
        self.infinum = sys.maxint

    def add_edge(self, node1, node2, seq_id):
        edge = Edge(node1, node2, seq_id)
        self.nodes[node1].edges.append(edge)
        self.nodes[node2].edges.append(edge)
        return edge

    def get_edges(self, node1, node2):
        edges = self.nodes[node1].edges
        return filter(lambda e: e.left_node == node2 or e.right_node == node2, edges)

    def get_black_edges(self, node1, node2):
        edges = self.get_edges(node1, node2)
        return filter(lambda e: e.seq_id is None, edges)

    def get_colored_edges(self, node1, node2):
        edges = self.get_edges(node1, node2)
        return filter(lambda e: e.seq_id is not None, edges)

    def remove_edges(self, node1, node2):
        to_del = self.get_edges(node1, node2)
        self.nodes[node1].edges = [e for e in self.nodes[node1].edges
                                        if e not in to_del]
        self.nodes[node2].edges = [e for e in self.nodes[node2].edges
                                        if e not in to_del]

    def neighbors(self, node):
        neighbors = set()
        for edge in self.nodes[node].edges:
            other_node = edge.left_node if edge.left_node != node else edge.right_node
            neighbors.add(other_node)

        return list(neighbors)

    def is_bifurcation(self, node):
        neighbors = self.neighbors(node)
        if len(neighbors) > 2:
            return True

        #all edges should be either black or colored
        for neighbor in neighbors:
            edges = self.get_edges(node, neighbor)
            seq_ids = map(lambda e: e.seq_id, edges)
            if None in seq_ids and len(seq_ids) > 1:
                return True

        return False

    def get_fragmented_blocks(self):
        set_objects = {}
        def get_set_obj(edge):
            if edge not in set_objects:
                set_objects[edge] = SetObj(edge)
            return set_objects[edge]

        for node in self.nodes:
            if not self.is_bifurcation(node):
                continue

            neighbors = self.neighbors(node)
            if len(neighbors) != 3:
                continue
            if self.infinum not in neighbors:
                continue
            neighbors.remove(self.infinum)

            left_node, right_node = neighbors
            left_black = filter(lambda e: e.seq_id is None,
                                self.nodes[left_node].edges)[0]
            right_black = filter(lambda e: e.seq_id is None,
                                 self.nodes[right_node].edges)[0]
            Union(get_set_obj(left_black), get_set_obj(right_black))

        groups_dict = defaultdict(list)
        for obj in set_objects.itervalues():
            groups_dict[Find(obj).obj].append(obj.obj)
        return groups_dict

    def get_permutations(self):
        ########
        next_edge_id = [1]
        edge_ids = {}
        def get_edge_id(edge):
            if not edge in edge_ids:
                edge_ids[edge] = next_edge_id[0]
                next_edge_id[0] += 1
            return edge_ids[edge]
        ########

        permutations = {}
        for seq_id, edge in self.origins.iteritems():
            blocks = []
            prev = edge
            edge = edge.next_edge

            while edge is not None:
                black_edges = self.get_black_edges(prev.right_node, edge.left_node)
                black_edge = black_edges[0]

                block_id = get_edge_id(black_edge)
                sign = 1 if black_edge.right_node == edge.left_node else -1
                start = prev.right_pos
                end = edge.left_pos
                length = end - start
                assert length >= 0

                blocks.append(Block(sign * block_id, start, length))
                prev, edge = edge, edge.next_edge

            permutations[seq_id] = blocks

        enum_groups = {}
        for repr_edge, edges in self.get_fragmented_blocks().iteritems():
            enum_groups[get_edge_id(repr_edge)] = map(get_edge_id, edges)

        return permutations, enum_groups


def build_graph(permutations):
    graph = BreakpointGraph()
    for seq_id, blocks in permutations.iteritems():
        #black edges
        for block in blocks:
            abs_block = abs(block.id)
            if not graph.get_black_edges(abs_block, -abs_block):
                graph.add_edge(abs_block, -abs_block, None)

        #chromosome ends
        head_edge = graph.add_edge(graph.infinum, blocks[0].id, seq_id)
        head_edge.right_pos = blocks[0].start
        tail_edge = graph.add_edge(-blocks[-1].id, graph.infinum, seq_id)
        tail_edge.left_pos = blocks[-1].start + blocks[-1].length
        graph.origins[seq_id] = head_edge

        prev_edge = head_edge
        edge = head_edge

        #adjacencies
        for block1, block2 in izip(blocks[:-1], blocks[1:]):
            left_pos = block1.start + block1.length
            right_pos = block2.start
            if right_pos < left_pos:
                print("WARNING: overlapping blocks")
                print(block1.start, block1.start + block1.length,
                      block2.start, block2.start + block2.length,
                      "|", right_pos - left_pos, seq_id)
            edge = graph.add_edge(-block1.id, block2.id, seq_id)
            edge.left_pos, edge.right_pos = left_pos, right_pos
            prev_edge.next_edge = edge
            edge.prev_edge = prev_edge
            prev_edge = edge

        edge.next_edge = tail_edge
        tail_edge.prev_edge = edge

    return graph


#################
class SetObj:
    def __init__(self, obj):
        self.obj = obj
        MakeSet(self)

def MakeSet(x):
     x.parent = x
     x.rank   = 0

def Union(x, y):
     xRoot = Find(x)
     yRoot = Find(y)
     if xRoot.rank > yRoot.rank:
         yRoot.parent = xRoot
     elif xRoot.rank < yRoot.rank:
         xRoot.parent = yRoot
     elif xRoot != yRoot:
         yRoot.parent = xRoot
         xRoot.rank = xRoot.rank + 1

def Find(x):
     if x.parent == x:
        return x
     else:
        x.parent = Find(x.parent)
        return x.parent
