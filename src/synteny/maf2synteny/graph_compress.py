import breakpoint_graph

#extend non-branching paths
def extend_path(graph, prev_node, cur_node, max_gap):
    path = [prev_node, cur_node]
    while True:
        #check if it's a bufurcation point
        if graph.is_bifurcation(cur_node) or graph.infinum in [prev_node, cur_node]:
            break

        #check distance
        edges = graph.get_colored_edges(prev_node, cur_node)
        get_len = lambda e: abs(e.right_pos - e.left_pos)
        gaps = map(get_len, edges)
        if gaps and max(gaps) > max_gap:
            break

        neighbors = graph.neighbors(cur_node)
        other_node = neighbors[0] if neighbors[0] != prev_node else neighbors[1]
        cur_node, prev_node = other_node, cur_node
        path.append(cur_node)

    return path


def compress_path(graph, path):
    if len(path) < 3:
        return False

    #ensure we start and end with black edge
    if not graph.get_black_edges(path[0], path[1]):
        del path[0]
    if not graph.get_black_edges(path[-2], path[-1]):
        del path[-1]

    if len(path) == 2:
        return False

    #updating graph
    graph.remove_edges(path[0], path[1])
    graph.remove_edges(path[-2], path[-1])
    graph.add_edge(path[0], path[-1], None)

    #updating links
    adjacencies = graph.get_colored_edges(path[1], path[2])
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


def compress_graph(graph, max_gap):
    #visited = set()
    num_compressed = 0
    to_delete = []

    for node in graph.nodes:
        if not graph.is_bifurcation(node):
            continue

        for neighbor in graph.neighbors(node):
            path = extend_path(graph, node, neighbor, max_gap)
            if compress_path(graph, path):
                num_compressed += 1
                to_delete.extend(path[1:-1])


    for node in to_delete:
        del graph.nodes[node]

    return num_compressed
