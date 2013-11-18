import networkx as nx
from permutation import *
from debug import DebugConfig, write_dot
import phylogeny as phylo


Connection = namedtuple("Connection", ["start", "end"])


class BreakpointGraph:
    def __init__(self):
        self.bp_graph = nx.MultiGraph()
        self.targets = []
        self.references = []


    def build_from(self, perm_container, circular):
        #TODO: add target here
        for perm in perm_container.ref_perms_filtered:
            prev = None
            for block in perm.iter_blocks(circular):
                if not prev:
                    prev = block
                    continue

                left_block = prev
                right_block = block

                self.bp_graph.add_node(-left_block)
                self.bp_graph.add_node(right_block)
                self.bp_graph.add_edge(-left_block, right_block, ref_id=perm.ref_id)

                prev = block

        for perm in perm_container.target_perms_filtered:
            if perm.ref_id not in self.targets:
                self.targets.append(perm.ref_id)
        for perm in perm_container.ref_perms_filtered:
            if perm.ref_id not in self.references:
                self.references.append(perm.ref_id)


    def find_adjacencies(self, phylogeny):
        adjacencies = {}
        subgraphs = nx.connected_component_subgraphs(self.bp_graph)
        for comp_id, subgraph in enumerate(subgraphs):
            #TODO: check for trivial case here
            weighted_graph = make_weighted(subgraph, phylogeny, self.targets, self.references)
            chosen_edges = split_graph(weighted_graph)

            for edge in chosen_edges:
                adjacencies[-edge[0]] = Connection(-edge[0], edge[1])
                adjacencies[-edge[1]] = Connection(-edge[1], edge[0])

            if DebugConfig.get_writer().debugging:
                debug_dir = DebugConfig.get_writer().debug_dir
                debug_draw_component(comp_id, weighted_graph, subgraph, debug_dir)

        return adjacencies


def debug_draw_component(comp_id, weighted_graph, breakpoint_graph, debug_dir):
    if len(breakpoint_graph) == 2:
        return

    for e in weighted_graph.edges_iter():
        weighted_graph[e[0]][e[1]]["label"] = ("{0:7.4f}"
                                    .format(weighted_graph[e[0]][e[1]]["weight"]))
    bg_out = os.path.join(debug_dir, "comp{0}-bg.dot".format(comp_id))
    weighted_out = os.path.join(debug_dir, "comp{0}-weighted.dot".format(comp_id))
    write_dot(breakpoint_graph, open(bg_out, "w"))
    write_dot(weighted_graph, open(weighted_out, "w"))


def split_graph(graph):
    #since we want minimum weight matching
    for e in graph.edges_iter():
        graph[e[0]][e[1]]["weight"] = -graph[e[0]][e[1]]["weight"]

    edges = nx.max_weight_matching(graph, maxcardinality=True)

    #nx algorithm returns duplicated adacency pairs
    unique_edges = []
    for v1, v2 in edges.iteritems():
        if not (v2, v1) in unique_edges:
            unique_edges.append((v1, v2))

    return unique_edges


def make_weighted(graph, phylogeny, target_ids, ref_ids):
    g = nx.Graph()
    g.add_nodes_from(graph.nodes())

    #trivial case
    if len(graph) == 2:
        node_1, node_2 = graph.nodes()
        g.add_edge(node_1, node_2, weight=1)
        return g

    #non-trivial
    for node in graph.nodes():
        adjacencies = {}
        for neighbor in graph.neighbors(node):
            for edge in graph[node][neighbor].values():
                adjacencies[edge["ref_id"]] = neighbor

        for ref_id in ref_ids:
            if ref_id not in adjacencies:
                adjacencies[ref_id] = None  #"void" state in paper

        target_id = target_ids[0]
        for neighbor in graph.neighbors(node):
            adjacencies[target_id] = neighbor
            breaks_weight, nbreaks, tree = phylogeny.estimate_tree(adjacencies)

            update_edge(g, node, neighbor, breaks_weight)

    return g


def update_edge(graph, v1, v2, weight):
    if not graph.has_edge(v1, v2):
        graph.add_edge(v1, v2, weight=weight)
    else:
        graph[v1][v2]["weight"] += weight

