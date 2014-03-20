from itertools import combinations
from collections import defaultdict

import breakpoint_graph
from graph_compress import extend_path

def collapse_bulge(graph, branches):
    """
    print map(str, branch1), map(str, branch2)

    for branch in [branch1, branch2]:
        print "0 -- 1", map(str, graph.get_edges(branch[0], branch[1]))
        if len(branch) > 3:
            print "1 -- 2", map(str, graph.get_edges(branch[1], branch[2]))
            print "3 -- 4", map(str, graph.get_edges(branch[2], branch[3]))
        print ""
    """
    for branch in branches:
        if len(branch) not in [2, 4]:
            return False

    for branch in branches:
        assert graph.infinum not in branch
        assert len(branch) in [2, 4]

        if len(branch) == 2:
            continue

        assert len(branch) == 4
        for adj in graph.get_colored_edges(branch[0], branch[1]):
            next_adj = adj.next_edge
            new_adj = None
            if next_adj.left_node in branch:
                next_adj = next_adj.next_edge
                prev_adj = adj.prev_edge

                new_adj = graph.add_edge(branch[0], branch[-1], adj.seq_id)
                new_adj.left_pos = adj.left_pos
                new_adj.right_pos = adj.next_edge.right_pos
            else:
                prev_adj = adj.prev_edge.prev_edge

                new_adj = graph.add_edge(branch[-1], branch[0], adj.seq_id)
                new_adj.left_pos = adj.prev_edge.left_pos
                new_adj.right_pos = adj.right_pos

            new_adj.next_edge = next_adj
            new_adj.prev_edge = prev_adj

            next_adj.prev_edge = new_adj
            prev_adj.next_edge = new_adj

        graph.remove_edges(branch[0], branch[1])
        graph.remove_edges(branch[-2], branch[-1])
    return True


def remove_bulges(graph, max_gap):
    num_collapsed = 0

    for node in graph.nodes:
        if not graph.is_bifurcation(node):
            continue

        by_end = defaultdict(list)
        for neighbor in graph.neighbors(node):
            path = extend_path(graph, node, neighbor, max_gap)
            by_end[path[-1]].append(path)

        #checking bulge structure: -<=>-
        if len(by_end) != 2:
            continue
        branches = None
        left_flank = None
        for b in by_end.values():
            if len(b) == 1:
                left_flank = graph.get_black_edges(b[0][0], b[0][1])
            else:
                branches = b
        if not (branches and left_flank):
            continue

        path_end = branches[0][-1]
        branch_repr = [b[-2] for b in branches]
        other_neighbors = []
        for node in graph.neighbors(path_end):
            if node not in branch_repr:
                other_neighbors.append(node)
        if len(other_neighbors) > 1 or not other_neighbors:
            continue
        right_flank = graph.get_black_edges(path_end, other_neighbors[0])
        if not right_flank:
            continue
        ##############

        if collapse_bulge(graph, branches):
            num_collapsed += 1

    return num_collapsed
