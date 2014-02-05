import sys
import networkx as nx
from collections import defaultdict
from itertools import izip
from Queue import Queue

from permutations import Block

def get_lcb(permutations, block_size):
    graph = build_graph(permutations)
    compress_graph(graph)
    return graph.get_permutations()


class Edge:
    def __init__(self, left_node, right_node, genome):
        self.genome = genome
        self.left_node = left_node
        self.right_node = right_node
        self.prev_edge = None
        self.next_edge = None


    def has_node(self, node):
        return self.left_node == node or self.right_node == node


    def __str__(self):
        left = str(self.left_node) if self.left_node != sys.maxint else "inf"
        right = str(self.right_node) if self.right_node != sys.maxint else "inf"
        return "({0}, {1}, {2})".format(left, right, self.genome)


class Node:
    def __init__(self):
        self.edges = set()
        self.coordinates = {}


class BreakpointGraph:
    def __init__(self):
        self.origins = {}
        self.nodes = defaultdict(Node)
        self.infinum = sys.maxint


    def add_edge(self, node1, node2, genome):
        edge = Edge(node1, node2, genome)
        self.nodes[node1].edges.add(edge)
        self.nodes[node2].edges.add(edge)
        return edge

    def has_black_edge(self, node1, node2):
        edges = self.get_edges(node1, node2)
        edges = filter(lambda e: e.genome == None, edges)
        return len(edges) > 0


    def get_edges(self, node1, node2):
        edges = self.nodes[node1].edges
        return filter(lambda e: e.left_node == node2 or e.right_node == node2, edges)


    def remove_edges(self, node1, node2):
        edges = self.get_edges(node1, node2)
        for e in edges:
            self.nodes[node1].edges.remove(e)
            self.nodes[node2].edges.remove(e)


    def neighbours(self, node):
        neighbours = set()
        for edge in self.nodes[node].edges:
            other_node = edge.left_node if edge.left_node != node else edge.right_node
            neighbours.add(other_node)

        return list(neighbours)


    #extend non-branching paths
    def extend_path(self, prev_node, cur_node):
        path = [prev_node, cur_node]
        while True:
            neighbours = self.neighbours(cur_node)
            if len(neighbours) != 2 or cur_node == self.infinum:
                break

            other_node = neighbours[0] if neighbours[0] != prev_node else neighbours[1]
            cur_node, prev_node = other_node, cur_node
            path.append(cur_node)

        return path


    def compress_path(self, path):
        if len(path) < 3:
            return

        #ensure we start and end with black edge
        if self.get_edges(path[0], path[1])[0].genome != None:
            del path[0]
        if self.get_edges(path[-2], path[-1])[0].genome != None:
            del path[-1]

        if len(path) == 2:
            return

        #print map(str, path)

        #updating graph
        self.remove_edges(path[0], path[1])
        self.remove_edges(path[-2], path[-1])
        self.add_edge(path[0], path[-1], None)

        #updating links
        adjacencies = self.get_edges(path[1], path[2])
        for adj in adjacencies:
            head_adj = adj.next_edge
            while not head_adj.has_node(path[0]) and not head_adj.has_node(path[-1]):
                head_adj = head_adj.next_edge

            tail_adj = adj.prev_edge
            while not tail_adj.has_node(path[0]) and not tail_adj.has_node(path[-1]):
                tail_adj = tail_adj.prev_edge

            assert head_adj.genome == tail_adj.genome
            head_adj.prev_edge = tail_adj
            tail_adj.next_edge = head_adj


    def get_permutations(self):
        next_id = [1]
        block_ids = {}
        def get_id(edge):
            if not edge in block_ids:
                block_ids[edge] = next_id[0]
                next_id[0] += 1
            return block_ids[edge]

        permutations = {}
        for genome, edge in self.origins.iteritems():
            blocks = []
            prev = edge
            edge = edge.next_edge

            while edge is not None:
                all_edges = self.get_edges(prev.right_node, edge.left_node)
                black_edge = filter(lambda e: e.genome is None, all_edges)[0]
                assert len(filter(lambda e: e.genome is None, all_edges)) == 1

                block_id = get_id(black_edge)
                sign = 1 if black_edge.right_node == edge.left_node else -1
                start = self.nodes[black_edge.left_node].coordinates[genome]
                end = self.nodes[black_edge.right_node].coordinates[genome]
                if start > end:
                    start, end = end, start
                length = end - start

                blocks.append(Block(sign * block_id, start, length))
                prev, edge = edge, edge.next_edge

            permutations[genome] = blocks

        return permutations


def build_graph(permutations):
    graph = BreakpointGraph()
    for genome, blocks in permutations.iteritems():
        #black edges
        for block in blocks:
            abs_block = abs(block.id)
            if not graph.has_black_edge(abs_block, -abs_block):
                graph.add_edge(abs_block, -abs_block, None)
            graph.nodes[abs_block].coordinates[genome] = block.start
            graph.nodes[-abs_block].coordinates[genome] = block.start + block.length
            #print block.start, block.length

        #chromosome ends
        head_edge = graph.add_edge(graph.infinum, blocks[0].id, genome)
        tail_edge = graph.add_edge(-blocks[-1].id, graph.infinum, genome)
        graph.origins[genome] = head_edge

        prev_edge = head_edge
        edge = head_edge

        #adjacencies
        for block1, block2 in izip(blocks[:-1], blocks[1:]):
            edge = graph.add_edge(-block1.id, block2.id, genome)
            prev_edge.next_edge = edge
            edge.prev_edge = prev_edge
            prev_edge = edge

        edge.next_edge = tail_edge
        tail_edge.prev_edge = edge

    return graph


def compress_graph(graph):
    opened_nodes = set([graph.infinum])
    closed_nodes = set()

    queue = Queue()
    queue.put(graph.infinum)

    while not queue.empty():
        node = queue.get()
        neighbours = graph.neighbours(node)

        for next_node in neighbours:
            path = graph.extend_path(node, next_node)
            path_end = path[-1]
            if not path_end in closed_nodes:
                graph.compress_path(path)

            if not path_end in opened_nodes:
                queue.put(path_end)
                opened_nodes.add(path_end)

        closed_nodes.add(node)
