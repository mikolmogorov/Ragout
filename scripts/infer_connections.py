from collections import namedtuple
from itertools import combinations, product
from breakpoint_graph import get_component_of, vertex_distance

Connection = namedtuple("Connection", ["start", "end", "distance"])



def case_trivial(graph, component, connected_comps, contig_index, num_ref):
    """
    a -- a
    """
    MIN_REF_THRESHOLD = 1   #TODO: think about it
    if len(component) != 2:
        return []

    num_edges = len(graph[component[0]].edges)
    if num_edges not in range(MIN_REF_THRESHOLD, num_ref + 1):
        return []

    #TODO: small contigs can be missed here!
    for fst, snd in [(0, 1), (1, 0)]:
        if abs(component[fst]) in contig_index and abs(component[snd]) not in contig_index:
            pair_comp = get_component_of(connected_comps, -component[snd])
            pair_id = pair_comp.index(-component[snd])
            other_id = abs(1 - pair_id)
            if pair_comp[other_id] in contig_index:

                print "indel found!"
                #TODO: distance estimation
                return [Connection(component[fst], pair_comp[other_id], 0)]

    if abs(component[0]) in contig_index and abs(component[1]) in contig_index:
        start = component[0]
        end = component[1]
        distance = vertex_distance(graph, start, end)
        return [Connection(start, end, distance)]

    return []


def all_equal(lst):
    if not lst:
        return True
    return lst.count(lst[0]) == len(lst)


def case_indel(graph, component, connected_comps, contig_index, num_ref):
    """
    a    -b
    |  \  |
    b     c
    """
    if len(component) != 4:
        return []

    found = False
    for v1, v2 in combinations(component, 2):
        if v1 == -v2:
            found = True
            similar = [v1, v2]
            different = filter(lambda v: v != v1 and v != v2, component)

    if not found:
        return []

    #checking topology
    #both b and -b have degree 1
    vertexes_1 = map(lambda e: e.vertex, graph[similar[0]].edges)
    vertexes_2 = map(lambda e: e.vertex, graph[similar[1]].edges)
    if not all_equal(vertexes_1) or not all_equal(vertexes_2):
        print "strange component!"
        return []
    #and they are connected to different vertexes
    if vertexes_1[0] == vertexes_2[0]:
        print "strange component!"
        return []

    if abs(similar[0]) in contig_index:
        #print "deletion in some references"
        connections = []
        for s in similar:
            assert abs(graph[s].edges[0].vertex) in contig_index
            distance = vertex_distance(graph, s, graph[s].edges[0].vertex)
            connections.append(Connection(s, graph[s].edges[0].vertex, distance))
        return connections
    else:
        #print "deletion in assembly and references"
        distance = vertex_distance(graph, different[0], different[1])
        assert abs(different[0]) in contig_index and abs(different[1]) in contig_index
        return [Connection(different[0], different[1], distance)]


def case_double_indel(graph, component, connected_comps, contig_index, num_ref):
    #TODO: check graph topology
    if len(component) != 6:
        return []

    comp_sorted = sorted(component, key=abs)
    pairs = []
    for i in xrange(len(comp_sorted) - 1):
        if abs(comp_sorted[i]) == abs(comp_sorted[i + 1]):
            pairs.append((comp_sorted[i], comp_sorted[i + 1]))
    if len(pairs) != 2:
        return []

    others = filter(lambda v: v not in pairs[0] + pairs[1], comp_sorted)
    print pairs, others

    if abs(others[0]) not in contig_index or abs(others[1]) not in contig_index:
        print "indel borders are not in contigs!"
        #assert False
        return []

    for i, j in [(0, 1), (1, 0)]:
        if abs(pairs[i][0]) in contig_index and abs(pairs[j][0]) not in contig_index:
            print "indel of two blocks, with one presented in assembly"
            edges1 = filter(lambda e: abs(e.vertex) == abs(pairs[i][0]), graph[others[0]].edges)
            v1 = edges1[0].vertex
            edges2 = filter(lambda e: abs(e.vertex) == abs(pairs[i][0]), graph[others[1]].edges)
            v2 = edges2[0].vertex
            dist1 = vertex_distance(graph, others[0], v1)
            dist2 = vertex_distance(graph, others[1], v2)
            return [Connection(others[0], v1, dist1), Connection(others[1], v2, dist2)]

    if abs(pairs[0][0]) in contig_index and abs(pairs[1][0]) in contig_index:
        print "indel of two blocks, with both presented in assembly"
        conns = []
        for i, j in product(pairs[0], pairs[1]):
            if j in map(lambda e: e.vertex, graph[i].edges):
                break

        dist = vertex_distance(graph, i, j)
        conns.append(Connection(i, j, dist))
        out = filter(lambda e: e.vertex in [-i, -j], graph[others[0]].edges)
        assert len(out) >= 1 and all_equal(map(lambda e: e.vertex, out))
        dist = vertex_distance(graph, others[0], out[0].vertex)
        conns.append(Connection(others[0], out[0].vertex, dist))

        out = filter(lambda e: e.vertex in [-i, -j], graph[others[1]].edges)
        assert len(out) >= 1 and all_equal(map(lambda e: e.vertex, out))
        dist = vertex_distance(graph, others[1], out[0].vertex)
        conns.append(Connection(others[1], out[0].vertex, dist))
        return conns

    if abs(pairs[0][0]) not in contig_index and abs(pairs[1][0]) not in contig_index:
        print "double deletion"
        return [Connection(others[0], others[1], 0)]

    return []


def simple_connections(graph, connected_comps, sibelia_output):
    contigs = sibelia_output.contigs
    contig_index = sibelia_output.build_contig_index()
    num_ref = sibelia_output.get_reference_count()
    connections = {}

    for component in connected_comps:
        conn_trivial = case_trivial(graph, component, connected_comps, contig_index, num_ref)
        conn_indel = case_indel(graph, component, connected_comps, contig_index, num_ref)
        conn_double_indel = case_double_indel(graph, component, connected_comps, contig_index, num_ref)

        for c in conn_trivial + conn_indel + conn_double_indel:
            assert abs(c.start) in contig_index
            assert abs(c.end) in contig_index
            connections[-c.start] = Connection(-c.start, c.end, c.distance)
            connections[-c.end] = Connection(-c.end, c.start, c.distance)

    print "connections infered:", len(connections) / 2
    return connections

