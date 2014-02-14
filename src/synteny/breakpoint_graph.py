import sys
import networkx as nx
from collections import defaultdict
from itertools import izip, combinations
from Queue import Queue

from permutations import Block

def get_lcb(permutations, max_gap):
    graph = build_graph(permutations)

    update = True
    while update:
        paths = graph.compress(max_gap)
        bulges = graph.find_bulges(max_gap)
        update = bool(paths + bulges)
        sys.stderr.write("Iteration: {0} paths, {1} bulges\n"
                                    .format(paths, bulges))
    return graph.get_permutations()


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
        self.edges = set()


class BreakpointGraph:
    def __init__(self):
        self.origins = {}
        self.nodes = defaultdict(Node)
        self.infinum = sys.maxint


    def add_edge(self, node1, node2, seq_id):
        edge = Edge(node1, node2, seq_id)
        self.nodes[node1].edges.add(edge)
        self.nodes[node2].edges.add(edge)
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
        edges = self.get_edges(node1, node2)
        for e in edges:
            self.nodes[node1].edges.remove(e)
            self.nodes[node2].edges.remove(e)


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


    #extend non-branching paths
    def extend_path(self, prev_node, cur_node, max_gap):
        path = [prev_node, cur_node]
        while True:
            #check if it's a bufurcation point
            if self.is_bifurcation(cur_node) or self.infinum in [prev_node, cur_node]:
                break

            #check distance
            edges = self.get_colored_edges(prev_node, cur_node)
            get_len = lambda e: abs(e.right_pos - e.left_pos)
            gaps = map(get_len, edges)
            if gaps and max(gaps) > max_gap:
                break

            neighbors = self.neighbors(cur_node)
            other_node = neighbors[0] if neighbors[0] != prev_node else neighbors[1]
            cur_node, prev_node = other_node, cur_node
            path.append(cur_node)

        return path


    def compress_path(self, path):
        if len(path) < 3:
            return False

        #ensure we start and end with black edge
        if not self.get_black_edges(path[0], path[1]):
            del path[0]
        if not self.get_black_edges(path[-2], path[-1]):
            del path[-1]

        if len(path) == 2:
            return False

        #updating graph
        self.remove_edges(path[0], path[1])
        self.remove_edges(path[-2], path[-1])
        self.add_edge(path[0], path[-1], None)

        #updating links
        adjacencies = self.get_colored_edges(path[1], path[2])
        for adj in adjacencies:
            head_adj = adj.next_edge
            while not head_adj.has_node(path[0]) and not head_adj.has_node(path[-1]):
                head_adj = head_adj.next_edge

            tail_adj = adj.prev_edge
            while not tail_adj.has_node(path[0]) and not tail_adj.has_node(path[-1]):
                tail_adj = tail_adj.prev_edge

            assert head_adj.seq_id == tail_adj.seq_id
            head_adj.prev_edge = tail_adj
            tail_adj.next_edge = head_adj

        return True


    def get_permutations(self):
        next_id = [1]
        block_ids = {}
        def get_id(edge):
            if not edge in block_ids:
                block_ids[edge] = next_id[0]
                next_id[0] += 1
            return block_ids[edge]

        permutations = {}
        for seq_id, edge in self.origins.iteritems():
            blocks = []
            prev = edge
            edge = edge.next_edge

            while edge is not None:
                black_edges = self.get_black_edges(prev.right_node, edge.left_node)
                black_edge = black_edges[0]

                block_id = get_id(black_edge)
                sign = 1 if black_edge.right_node == edge.left_node else -1
                start = prev.right_pos
                end = edge.left_pos
                length = end - start
                assert length >= 0

                blocks.append(Block(sign * block_id, start, length))
                prev, edge = edge, edge.next_edge

            permutations[seq_id] = blocks

        return permutations


    def compress(self, max_gap):
        #visited = set()
        num_compressed = 0
        to_delete = []

        for node in self.nodes:
            if not self.is_bifurcation(node):
                continue

            for neighbor in self.neighbors(node):
                path = self.extend_path(node, neighbor, max_gap)
                if self.compress_path(path):
                    num_compressed += 1
                    to_delete.extend(path[1:-1])


        for node in to_delete:
            del self.nodes[node]

        return num_compressed


    def collapse_bulge(self, branch1, branch2):
        #TODO: think about this
        if (len(self.neighbors(branch1[0])) > 3 or
            len(self.neighbors(branch1[-1])) > 3):
            return False

        for branch in [branch1, branch2]:
            if self.infinum in branch:
                return False
            if len(branch) not in [2, 4]:
                return False
            if len(branch) == 2 and self.get_black_edges(branch[0], branch[1]):
                return False

        #print map(str, branch1), map(str, branch2)

        for branch in [branch1, branch2]:
            if len(branch) == 2:
                continue

            assert len(branch) == 4
            for adj in self.get_colored_edges(branch[0], branch[1]):
                next_adj = adj.next_edge
                new_adj = None
                if next_adj.left_node in branch:
                    next_adj = next_adj.next_edge
                    prev_adj = adj.prev_edge

                    new_adj = self.add_edge(branch[0], branch[-1], adj.seq_id)
                    new_adj.left_pos = adj.left_pos
                    new_adj.right_pos = adj.next_edge.right_pos
                else:
                    prev_adj = adj.prev_edge.prev_edge

                    new_adj = self.add_edge(branch[-1], branch[0], adj.seq_id)
                    new_adj.left_pos = adj.prev_edge.left_pos
                    new_adj.right_pos = adj.right_pos

                new_adj.next_edge = next_adj
                new_adj.prev_edge = prev_adj

                next_adj.prev_edge = new_adj
                prev_adj.next_edge = new_adj

            self.remove_edges(branch[0], branch[1])
            self.remove_edges(branch[-2], branch[-1])
        return True


    def find_bulges(self, max_gap):
        num_collapsed = 0

        for node in self.nodes:
            if not self.is_bifurcation(node):
                continue

            update = True
            while update:
                update = False
                paths = {}

                for neighbor in self.neighbors(node):
                    path = self.extend_path(node, neighbor, max_gap)
                    paths[neighbor] = path

                for path1, path2 in combinations(paths.itervalues(), 2):
                    if path1[-1] != path2[-1]:
                        continue


                    if self.collapse_bulge(path1, path2):
                        num_collapsed += 1
                        update = True
                        break

        return num_collapsed


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
            edge = graph.add_edge(-block1.id, block2.id, seq_id)
            edge.left_pos = block1.start + block1.length
            edge.right_pos = block2.start
            assert edge.right_pos >= edge.left_pos

            prev_edge.next_edge = edge
            edge.prev_edge = prev_edge
            prev_edge = edge

        edge.next_edge = tail_edge
        tail_edge.prev_edge = edge

    return graph
