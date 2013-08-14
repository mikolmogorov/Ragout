from collections import namedtuple
from itertools import combinations, product
from breakpoint_graph import get_component_of, vertex_distance, get_connected_components

Connection = namedtuple("Connection", ["start", "end", "distance"])
Pair = namedtuple("Pair", ["fst", "snd"])

def all_equal(lst):
        if not lst:
            return True
        return lst.count(lst[0]) == len(lst)


class AdjacencyFinder:
    def __init__(self, graph, sibelia_output):
        self.graph = graph
        self.sibelia_output = sibelia_output


    def in_assembly(self, block):
        return abs(block) in self.contig_index


    def case_trivial(self, component):
        """
        a -- a
        """
        #TODO: think about it
        MIN_REF_THRESHOLD = 1
        if len(component) != 2:
            return []

        vertex_pair = Pair(component[0], component[1])
        num_edges = len(self.graph[vertex_pair.fst].edges)
        if num_edges < MIN_REF_THRESHOLD:
            return []

        for fst, snd in [(0, 1), (1, 0)]:
            if self.in_assembly(vertex_pair[fst]) and not self.in_assembly(vertex_pair[snd]):
                pair_comp = get_component_of(self.connected_comps, -vertex_pair[snd])
                pair_id = pair_comp.index(-vertex_pair[snd])
                other_id = abs(1 - pair_id)
                if self.in_assembly(pair_comp[other_id]):

                    print "indel found!"
                    #TODO: distance estimation
                    return [Connection(vertex_pair[fst], pair_comp[other_id], None)]

        if self.in_assembly(vertex_pair.fst) and self.in_assembly(vertex_pair.snd):
            start = vertex_pair.fst
            end = vertex_pair.snd
            distance = vertex_distance(self.graph, start, end)
            return [Connection(start, end, distance)]

        return []


    def case_indel(self, component):
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
                similar = Pair(v1, v2)
                different = Pair(*filter(lambda v: v != v1 and v != v2, component))

        if not found:
            return []

        #checking topology
        #both b and -b have degree 1
        vertexes_1 = map(lambda e: e.vertex, self.graph[similar.fst].edges)
        vertexes_2 = map(lambda e: e.vertex, self.graph[similar.snd].edges)
        if not all_equal(vertexes_1) or not all_equal(vertexes_2):
            print "strange component!"
            return []
        #and they are connected to different vertexes
        if vertexes_1[0] == vertexes_2[0]:
            print "strange component!"
            return []

        if self.in_assembly(similar.fst):
            print "deletion in some references"
            connections = []
            for s in similar:
                pair_vertex = self.graph[s].edges[0].vertex
                assert self.in_assembly(pair_vertex)
                distance = vertex_distance(self.graph, s, pair_vertex)
                connections.append(Connection(s, pair_vertex, distance))
            return connections
        else:
            print "deletion in assembly and references"
            distance = vertex_distance(self.graph, different.fst, different.snd)
            assert self.in_assembly(different.fst) and self.in_assembly(different.snd)
            return [Connection(different.fst, different.snd, distance)]


    def case_double_indel(self, component):
        #TODO: check graph topology
        if len(component) != 6:
            return []

        comp_sorted = sorted(component, key=abs)
        pairs = []
        for i in xrange(len(comp_sorted) - 1):
            if abs(comp_sorted[i]) == abs(comp_sorted[i + 1]):
                pairs.append(Pair(comp_sorted[i], comp_sorted[i + 1]))
        if len(pairs) != 2:
            return []

        others = Pair(*filter(lambda v: v not in pairs[0] + pairs[1], comp_sorted))
        #print pairs, others

        if not self.in_assembly(others.fst) or not self.in_assembly(others.snd):
            print "indel borders are not in contigs!"
            #assert False
            return []

        for i, j in [(0, 1), (1, 0)]:
            if self.in_assembly(pairs[i].fst) and not self.in_assembly(pairs[j].fst):
                print "indel of two blocks, with one presented in assembly"
                edges1 = filter(lambda e: abs(e.vertex) == abs(pairs[i].fst), self.graph[others.fst].edges)
                edges2 = filter(lambda e: abs(e.vertex) == abs(pairs[i].fst), self.graph[others.snd].edges)
                assert len(edges1) != 0 or len(edges2) != 0
                if len(edges1) > 0:
                    v1 = edges1[0].vertex
                    dist1 = vertex_distance(self.graph, others.fst, v1)
                else:
                    v1 = -edges2[0].vertex
                    dist1 = None
                if len(edges2) > 0:
                    v2 = edges2[0].vertex
                    dist2 = vertex_distance(self.graph, others.snd, v2)
                else:
                    v2 = -edges1[0].vertex
                    dist2 = None

                return [Connection(others.fst, v1, dist1), Connection(others.snd, v2, dist2)]

        if self.in_assembly(pairs[0].fst) and self.in_assembly(pairs[1].fst):
            print "indel of two blocks, with both presented in assembly"
            conns = []
            for i, j in product(pairs[0], pairs[1]):
                if j in map(lambda e: e.vertex, self.graph[i].edges):
                    break

            dist = vertex_distance(self.graph, i, j)
            conns.append(Connection(i, j, dist))
            out = filter(lambda e: e.vertex in [-i, -j], self.graph[others.fst].edges)
            assert len(out) >= 1 and all_equal(map(lambda e: e.vertex, out))
            dist = vertex_distance(self.graph, others.fst, out[0].vertex)
            conns.append(Connection(others.fst, out[0].vertex, dist))

            out = filter(lambda e: e.vertex in [-i, -j], self.graph[others.snd].edges)
            assert len(out) >= 1 and all_equal(map(lambda e: e.vertex, out))
            dist = vertex_distance(self.graph, others.snd, out[0].vertex)
            conns.append(Connection(others.snd, out[0].vertex, dist))
            return conns

        if not self.in_assembly(pairs[0].fst) and not self.in_assembly(pairs[1].fst):
            print "double deletion"
            return [Connection(others.fst, others.snd, None)]

        return []


    def find_adjacencies(self):
        self.contig_index = self.sibelia_output.build_contig_index()
        self.num_ref = self.sibelia_output.get_reference_count()
        self.connected_comps = get_connected_components(self.graph)
        connections = {}

        for component in self.connected_comps:
            conn_trivial = self.case_trivial(component)
            conn_indel = self.case_indel(component)
            conn_double_indel = self.case_double_indel(component)

            for c in conn_trivial + conn_indel + conn_double_indel:
                assert self.in_assembly(c.start)
                assert self.in_assembly(c.end)
                connections[-c.start] = Connection(-c.start, c.end, c.distance)
                connections[-c.end] = Connection(-c.end, c.start, c.distance)

        print "connections infered:", len(connections) / 2
        #for f, to in connections.iteritems():
        #    print f, to
        return connections

