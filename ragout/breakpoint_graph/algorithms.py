#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
Various algorithms for breakpoint graph processing
"""

from collections import deque

import networkx as nx

def min_weight_matching(graph):
    """
    Finds a perfect matching with minimum weight
    """
    for v1, v2 in graph.edges_iter():
        graph[v1][v2]["weight"] = -graph[v1][v2]["weight"] #want minimum weght

    MIN_LOG_SIZE = 20
    if len(graph) > MIN_LOG_SIZE:
        logger.debug("Finding perfect matching for a component of "
                     "size {0}".format(len(graph)))
    edges = nx.max_weight_matching(graph, maxcardinality=True)
    unique_edges = set()
    for v1, v2 in edges.items():
        if not (v2, v1) in unique_edges:
            unique_edges.add((v1, v2))

    return list(unique_edges)


def alternating_cycle(graph, node_1, node_2, target_id):
    """
    Determines if there is a cycle of alternating colors
    that goes through the given edge
    """
    def get_genome_ids((u, v)):
        return list(map(lambda e: e["genome_id"], graph[u][v].values()))

    simple_graph = nx.Graph()
    for (u, v) in graph.edges_iter():
        simple_graph.add_edge(u, v)
    assert not simple_graph.has_edge(node_1, node_2)

    good_path = False
    for path in nx.all_simple_paths(simple_graph, node_1, node_2):
        #2-break or 3-break
        if len(path) % 2 == 1 or len(path) / 2 > 3:
            continue

        edges = list(zip(path[:-1], path[1:]))
        odd_colors = list(map(get_genome_ids, edges[0::2]))
        even_colors = list(map(get_genome_ids, edges[1::2]))

        if not all(map(lambda e: set(e) == set([target_id]), even_colors)):
            continue

        common_genomes = set(odd_colors[0])
        for edge_colors in odd_colors:
            common_genomes = common_genomes.intersection(edge_colors)

        if common_genomes:
            good_path = True
            break

    return good_path


def mark_preferred_edges(graph, trusted_adjacencies):
    """
    Try to bring edges that are supported by previous iteration
    to the front
    """
    def bfs_search(start_node, goal_node):
        queue = deque()
        visited = set([start_node])
        queue.append((start_node, 0))
        while len(queue):
            cur_node, cur_depth = queue.popleft()
            if cur_node == goal_node:
                return cur_depth

            for other_node in graph.neighbors(cur_node):
                if -other_node not in visited:
                    queue.append((-other_node, cur_depth + 1))
                    visited.add(-other_node)

    for node in graph.nodes():
        if len(graph.neighbors(node)) == 1 or node not in trusted_adjacencies:
            continue

        for neighbor in graph.neighbors(node):
            num_steps = bfs_search(neighbor, trusted_adjacencies[node])
            print(node, neighbor, num_steps)
