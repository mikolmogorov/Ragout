import networkx as nx
from permutation import *
import phylogeny as phylo


Connection = namedtuple("Connection", ["start", "end"])


class BreakpointGraph:
    def __init__(self):
        self.bp_graph = nx.MultiGraph()


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

        nx.draw_graphviz(self.bp_graph)
        #nx.write_dot(self.bp_graph, "g.dot")


    def find_adjacencies(self, phylogeny):
        adjacencies = {}
        #component_counter = 0
        subgraphs = nx.connected_component_subgraphs(self.bp_graph)
        for subgraph in subgraphs:
            #TODO: check for trivial case here
            weighted_graph = make_weighted(subgraph, phylogeny)

            chosen_edges = split_graph(weighted_graph)

            for edge in chosen_edges:
                adjacencies[-edge[0]] = Connection(-edge[0], edge[1])
                adjacencies[-edge[1]] = Connection(-edge[1], edge[0])

            #if self.debug_dir and len(subgraph) > 2:
            #    for e in graph.edges_iter():
            #        weighted_graph[e[0]][e[1]]["label"] = "{0:5.2f}".format(weighted_graph[e[0]][e[1]]["weight"])
            #    bg_out = open(os.path.join(self.debug_dir, "comp{0}-bg.dot".format(component_counter)), "w")
            #    prelinks_out = os.path.join(self.debug_dir, "comp{0}-prelinks.dot".format(component_counter))
            #    write_colored_dot(subgraph, bg_out)
            #    nx.draw_graphviz(weighted_graph)
            #    nx.write_dot(weighted_graph, prelinks_out)

            #component_counter += 1

        return adjacencies


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


def make_weighted(graph, phylogeny):
    g = nx.Graph()
    g.add_nodes_from(graph.nodes())

    #trivial case
    if len(graph) == 2:
        node_1, node_2 = graph.nodes()
        g.add_edge(node_1, node_2, weight=1)
        #g.add_edge(node_2, node_1, weight=1)
        return g

    #non-trivial
    for node in graph.nodes():
        adjacencies = {}
        for neighbor in graph.neighbors(node):
            for edge in graph[node][neighbor].values():
                adjacencies[edge["ref_id"]] = neighbor

        #max_likelihood = float("-inf")
        #max_tree = None

        for neighbor in graph.neighbors(node):
            adjacencies["target"] = neighbor
            breaks_weight, nbreaks, tree = phylogeny.estimate_tree(adjacencies)

            update_edge(g, node, neighbor, breaks_weight)

            #if likelihood > max_likelihood:
            #    max_likelihood = likelihood
            #    max_tree = tree

        #if self.debug_dir:
        #    debug_pref = os.path.join(self.debug_dir, "comp{0}-".format(component_counter))
        #    debug_file = open(debug_pref + "node_{0}.dot".format(node), "w")
        #    phylo.tree_to_dot(max_tree, debug_file)

    return g

def update_edge(graph, v1, v2, weight):
    if not graph.has_edge(v1, v2):
        graph.add_edge(v1, v2, weight=weight)
    else:
        graph[v1][v2]["weight"] += weight

    #TODO: graph output
    #def write_colored_dot(self, dot_file):
    #    Colors = ["red", "green", "blue", "yellow", "black"]
    #    def output_subgraph(subgraph):
    #        for edge in subgraph.edges(data=True):
    #            color = Colors[edge[2]["color"] - 1]
    #            dot_file.write("""{0} -- {1} [color = "{2}"];\n"""
    #                            .format(edge[0], edge[1], color))

    #    dot_file.write("graph {\n")
    #    output_subgraph(graph)
    #    dot_file.write("}\n")
