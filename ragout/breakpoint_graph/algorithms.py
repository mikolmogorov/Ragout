#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
Various algorithms for breakpoint graph processing
"""

from __future__ import print_function
from collections import deque

import networkx as nx

def min_weight_matching(graph, preferred_edges):
    """
    Finds a perfect matching with minimum weight
    """
    max_weight, min_weight = float("-inf"), float("inf")
    for v1, v2 in graph.edges_iter():
        graph[v1][v2]["weight"] = -graph[v1][v2]["weight"] #want minimum weght
        max_weight = max(max_weight, graph[v1][v2]["weight"])
        min_weight = min(min_weight, graph[v1][v2]["weight"])

    boost = max_weight - min_weight
    for v1, v2 in graph.edges_iter():
        if preferred_edges.get(v1, None) == v2:
            #print(boost)
            graph[v1][v2]["weight"] += boost + 1

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


def get_preferred_edges(graph, trusted_adjacencies, restricted_nodes):
    """
    Try to bring edges that are supported by previous iteration
    to the front
    """
    preferred_edges = {}

    def bfs_search(start_node, goal_node):
        #print("Search from", start_node, "to", goal_node)
        queue = deque()
        visited = set([start_node, -start_node])
        queue.append((-start_node, 0))
        while len(queue):
            cur_node, cur_depth = queue.popleft()
            #print("Cur", cur_node)
            if cur_node == goal_node:
                return cur_depth
            if abs(cur_node) in restricted_nodes:
                continue

            for other_node in graph.neighbors(cur_node):
                if -other_node not in visited:
                    #print("Other", other_node, "depth", cur_depth)
                    queue.append((-other_node, cur_depth + 1))
                    visited.add(-other_node)
                    visited.add(other_node)

        return float("inf")

    for node in graph.nodes():
        if len(graph.neighbors(node)) == 1 or node not in trusted_adjacencies:
            continue

        #print("Trusted adj: {0} -- {1}".format(node, trusted_adjacencies[node]))
        steps = {}
        for neighbor in graph.neighbors(node):
            num_steps = bfs_search(neighbor, trusted_adjacencies[node])
            steps[neighbor] = num_steps
        #print(steps)

        if sorted(steps.values())[0] < sorted(steps.values())[1]:
            preferred_edges[node] = neighbor
            preferred_edges[neighbor] = node

    return preferred_edges
