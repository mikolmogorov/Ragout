from collections import namedtuple

Colors = ["red", "green", "blue", "yellow", "black"]

def output_graph(graph, dot_file):
    dot_file.write("graph {\n")
    used_vertexes = set()
    for node_id, node in graph.iteritems():
        for edge in node.edges:
            if edge.vertex not in used_vertexes:
                dot_file.write("""{0} -- {1} [color = "{2}"];\n"""
                                .format(node_id, edge.vertex, Colors[edge.color]))
        used_vertexes.add(node_id)

    #for i in xrange(max(graph.keys()) + 1):
    #    dot_file.write("""{0} -- {1} [color = "black"];\n""".format(i, -i))
    dot_file.write("}")


def get_connected_components(graph):
    con_comp = []

    visited = set()
    def dfs(vertex, component):
        visited.add(vertex)
        component.append(vertex)
        for edge in graph[vertex].edges:
            if edge.vertex not in visited:
                dfs(edge.vertex, component)

    for vertex in graph:
        if vertex not in visited:
            con_comp.append([])
            dfs(vertex, con_comp[-1])

    return con_comp
