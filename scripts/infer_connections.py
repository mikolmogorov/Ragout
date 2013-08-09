import numpy
from collections import namedtuple
from itertools import combinations

Connection = namedtuple("Connection", ["start", "end", "distance"])

def median(vals_list):
    return numpy.median(vals_list)

def get_component_of(connected_comps, vertex):
    for con in connected_comps:
        if vertex in con:
            return con
    return None


def case_trivial(graph, component, connected_comps, contig_index, num_ref):
    """
    a -- a
    """
    MIN_REF_THRESHOLD = 2   #TODO: think about it
    if len(component) != 2:
        return []

    num_edges = len(graph[component[0]].edges)
    if num_edges not in range(MIN_REF_THRESHOLD, num_ref + 1):
        return []

    for fst, snd in [(0, 1), (1, 0)]:
        #TODO: test this case (never happened yet)
        if abs(component[fst]) in contig_index and abs(component[snd]) not in contig_index:
            pair_comp = get_component_of(connected_comps, -component[snd])
            pair_id = pair_comp.index(-component[snd])
            other_id = abs(1 - pair_id)
            if pair_comp[other_id] in contigs:
                print "indel found!"
                return [Connection(component[fst], pair_comp[other_id], None)]

    if abs(component[0]) in contig_index and abs(component[1]) in contig_index:
        start = component[0]
        end = component[1]
        distance = median(map(lambda e:e.distance, graph[start].edges))
        return [Connection(start, end, distance)]

    return []


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
    #TODO: check graph structure

    if abs(similar[0]) in contig_index:
        print "deletion in some references"
        connections = []
        for s in similar:
            distance = median(map(lambda e : e.distance, graph[s].edges))
            #if abs(s) in contig_index and abs(graph[s].edges[0].vertex) in contig_index:
            connections.append(Connection(s, graph[s].edges[0].vertex, distance))
        return connections
    else:
        print "deletion in assembly and (possibly) references"
        edges = filter(lambda e : e.vertex == different[1], graph[different[0]].edges)
        distance = median(map(lambda e : e.distance, edges))
        #if abs(different[0]) in contig_index and abs(different[1]) in contig_index:
        return [Connection(different[0], different[1], distance)]


def simple_connections(graph, connected_comps, sibelia_output):
    contigs = sibelia_output.contigs
    contig_index = sibelia_output.build_contig_index()
    num_ref = sibelia_output.get_reference_count()
    connections = {}

    for component in connected_comps:
        conn_trivial = case_trivial(graph, component, connected_comps, contig_index, num_ref)
        conn_indel = case_indel(graph, component, connected_comps, contig_index, num_ref)
        for c in conn_trivial + conn_indel:
            assert abs(c.start) in contig_index
            assert abs(c.end) in contig_index
            connections[-c.start] = Connection(-c.start, c.end, c.distance)
            connections[-c.end] = Connection(-c.end, c.start, c.distance)

    print "connections infered:", len(connections)
    return connections

