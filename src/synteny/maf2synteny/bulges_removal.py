from itertools import combinations

import breakpoint_graph
from graph_compress import extend_path

def collapse_bulge(graph, branch1, branch2):
    #TODO: think about this
    if (len(graph.neighbors(branch1[0])) > 3 or
        len(graph.neighbors(branch1[-1])) > 3):
        return False

    for branch in [branch1, branch2]:
        if graph.infinum in branch:
            return False
        if len(branch) not in [2, 4]:
            return False
        if len(branch) == 2 and graph.get_black_edges(branch[0], branch[1]):
            return False

    #print map(str, branch1), map(str, branch2)

    for branch in [branch1, branch2]:
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

        update = True
        while update:
            update = False
            paths = {}

            for neighbor in graph.neighbors(node):
                path = extend_path(graph, node, neighbor, max_gap)
                paths[neighbor] = path

            for path1, path2 in combinations(paths.itervalues(), 2):
                if path1[-1] != path2[-1]:
                    continue


                if collapse_bulge(graph, path1, path2):
                    num_collapsed += 1
                    update = True
                    break

    return num_collapsed
