from collections import namedtuple
from itertools import combinations, product
from breakpoint_graph import get_component_of, vertex_distance, get_connected_components

Connection = namedtuple("Connection", ["start", "end", "distance"])
Pair = namedtuple("Pair", ["fst", "snd"])

def all_equal(lst):
    if not lst:
        return True
    return lst.count(lst[0]) == len(lst)


"""
Deeply magical part of the tool
"""
class AdjacencyFinder:
    def __init__(self, graph, sibelia_output):
        self.graph = graph
        self.sibelia_output = sibelia_output


    def in_assembly(self, block):
        return self.graph[block].in_assembly


    def component_simple(self, components):
        """
        Cases with components of length 2
        """
        MIN_REF_THRESHOLD = 1

        connections = []
        resolved_comps_ids = []

        for comp_num, component in enumerate(components):
            if len(component) != 2:
                continue
            if comp_num in resolved_comps_ids:
                continue

            vertex_pair = Pair(*component)
            num_edges = len(self.graph[vertex_pair.fst].edges)
            if num_edges < MIN_REF_THRESHOLD:
                continue

            #one block from component is missing from assembly
            #could be indel
            ###############
            for fst, snd in [(0, 1), (1, 0)]:
                if self.in_assembly(vertex_pair[fst]) and not self.in_assembly(vertex_pair[snd]):
                    pair_comp = get_component_of(self.connected_comps, -vertex_pair[snd])
                    if len(pair_comp) != 2:
                        continue

                    pair_id = pair_comp.index(-vertex_pair[snd])
                    other_id = abs(1 - pair_id)

                    if self.in_assembly(pair_comp[other_id]):
                        print "indel found!"
                        start = vertex_pair[fst]
                        middle = vertex_pair[snd]
                        end = pair_comp[other_id]

                        distance = (vertex_distance(self.graph, start, middle) +
                                    vertex_distance(self.graph, -middle, end))
                        connections.append(Connection(start, end, distance))
                        resolved_comps_ids.extend([comp_num, components.index(pair_comp)])
                        continue
            #######

            #both exists in assembly
            ##############
            if self.in_assembly(vertex_pair.fst) and self.in_assembly(vertex_pair.snd):
                start = vertex_pair.fst
                end = vertex_pair.snd

                distance = vertex_distance(self.graph, start, end)
                connections.append(Connection(start, end, distance))
                resolved_comps_ids.append(comp_num)
            ######

        return (connections, resolved_comps_ids)


    def component_indel(self, components):
        """
        special case with 4 vertexes
        a    -b
        |  \  |
        b     c
        """
        connections = []
        resolved_comps_ids = []

        for comp_num, component in enumerate(components):
            if len(component) != 4:
                continue

            #checking topology
            ###################
            #should have two similar vertexes
            found = False
            for v1, v2 in combinations(component, 2):
                if v1 == -v2:
                    found = True
                    similar = Pair(v1, v2)
                    different = Pair(*filter(lambda v: v != v1 and v != v2, component))
            if not found:
                continue

            #both b and -b have degree 1
            vertexes_1 = map(lambda e: e.vertex, self.graph[similar.fst].edges)
            vertexes_2 = map(lambda e: e.vertex, self.graph[similar.snd].edges)
            if not all_equal(vertexes_1) or not all_equal(vertexes_2):
                print "strange component!"
                continue
            #and they are connected to different vertexes
            if vertexes_1[0] == vertexes_2[0]:
                print "strange component!"
                continue
            ####################
            if not self.in_assembly(different.fst) or not self.in_assembly(different.snd):
                print "indel borders are not in assembly"
                continue

            #checking, if indel has occured in assembly
            ###########
            if self.in_assembly(similar.fst):
                print "indel in some references"
                for s in similar:
                    pair_vertex = self.graph[s].edges[0].vertex

                    distance = vertex_distance(self.graph, s, pair_vertex)
                    connections.append(Connection(s, pair_vertex, distance))
                resolved_comps_ids.append(comp_num)
            else:
                print "deletion in assembly and references"
                distance = vertex_distance(self.graph, different.fst, different.snd)
                connections.append(Connection(different.fst, different.snd, distance))
                resolved_comps_ids.append(comp_num)
            ############

        return (connections, resolved_comps_ids)


    #TODO: check graph topology
    def component_double_indel(self, components):
        """
        case with component of size 6 (double indel)
        """
        connections = []
        resolved_comps_ids = []

        for comp_num, component in enumerate(components):
            if len(component) != 6:
                continue

            #findig pairs of oppisite vertexes (a and -a)
            ###############
            comp_sorted = sorted(component, key=abs)
            pairs = []
            for i in xrange(len(comp_sorted) - 1):
                if abs(comp_sorted[i]) == abs(comp_sorted[i + 1]):
                    pairs.append(Pair(comp_sorted[i], comp_sorted[i + 1]))
            if len(pairs) != 2:
                continue

            others = Pair(*filter(lambda v: v not in pairs[0] + pairs[1], comp_sorted))
            ############

            #check if "unique" blocks are presented in assembly
            if not self.in_assembly(others.fst) or not self.in_assembly(others.snd):
                print "indel borders missing from contigs!"
                continue

            #case with one of two blocks presented in assembly
            ################
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

                    connections.extend([Connection(others.fst, v1, dist1),
                                        Connection(others.snd, v2, dist2)])
                    resolved_comps_ids.append(comp_num)
            ##################

            #each of two blocks is presented in assembly
            #############
            if self.in_assembly(pairs[0].fst) and self.in_assembly(pairs[1].fst):
                print "indel of two blocks, with both presented in assembly"
                for i, j in product(pairs[0], pairs[1]):
                    if j in map(lambda e: e.vertex, self.graph[i].edges):
                        break

                dist = vertex_distance(self.graph, i, j)
                connections.append(Connection(i, j, dist))

                out = filter(lambda e: e.vertex in [-i, -j], self.graph[others.fst].edges)
                assert len(out) >= 1 and all_equal(map(lambda e: e.vertex, out))
                dist = vertex_distance(self.graph, others.fst, out[0].vertex)
                connections.append(Connection(others.fst, out[0].vertex, dist))

                out = filter(lambda e: e.vertex in [-i, -j], self.graph[others.snd].edges)
                assert len(out) >= 1 and all_equal(map(lambda e: e.vertex, out))
                dist = vertex_distance(self.graph, others.snd, out[0].vertex)
                connections.append(Connection(others.snd, out[0].vertex, dist))

                resolved_comps_ids.append(comp_num)
                continue
            ################

            #each of two blocks is missing
            if not self.in_assembly(pairs[0].fst) and not self.in_assembly(pairs[1].fst):
                print "double deletion"
                connections.append(Connection(others.fst, others.snd, None))
                resolved_comps_ids.append(comp_num)
            #################

        return (connections, resolved_comps_ids)

    def component_substitution(self, components):
        """
        substitution of one block, corresponds to two components of length 3
        """
        connections = []
        resolved_comps_ids = []

        for comp_num, component in enumerate(components):
            if len(component) != 3:
                continue
            if comp_num in resolved_comps_ids:
                continue

            #finding vertexes with degree 1
            variant_vertexes = []
            for vertex in component:
                if all_equal(map(lambda e: e.vertex, self.graph[vertex].edges)):
                    variant_vertexes.append(vertex)
            if len(variant_vertexes) != 2:
                continue
            common_vertex = filter(lambda v: v not in variant_vertexes, component)[0]

            #find corresponding component
            corr_comp = get_component_of(components, -variant_vertexes[0])
            if not corr_comp or len(corr_comp) != 3:
                continue

            #check component topology
            for vertex in variant_vertexes:
                if not all_equal(map(lambda e: e.vertex, self.graph[-vertex].edges)):
                    continue
            ############

            #determine, which of variant blocks is presented in essembly
            for i in [0, 1]:
                if self.in_assembly(variant_vertexes[i]):
                    break
            if self.in_assembly(variant_vertexes[abs(1 - i)]):
                print "strange substitution!"
                continue
            ##########
            print "substitution found"
            start = common_vertex
            middle = variant_vertexes[i]
            end = self.graph[-variant_vertexes[i]].edges[0].vertex

            dist1 = vertex_distance(self.graph, start, middle)
            connections.append(Connection(start, middle, dist1))
            dist2 = vertex_distance(self.graph, -middle, end)
            connections.append(Connection(-middle, end, dist2))
            resolved_comps_ids.extend([comp_num, components.index(corr_comp)])

        return (connections, resolved_comps_ids)


    def find_adjacencies(self):
        self.num_ref = self.sibelia_output.get_reference_count()
        self.connected_comps = get_connected_components(self.graph)
        connections = {}

        functions = [AdjacencyFinder.component_simple,
                     AdjacencyFinder.component_indel,
                     AdjacencyFinder.component_double_indel,
                     AdjacencyFinder.component_substitution]
        components = self.connected_comps
        result = []

        #print components
        for fun in functions:
            connect, resolved_comps_ids = fun(self, components)
            result.extend(connect)
            filtered = filter(lambda (i, _): i not in resolved_comps_ids, enumerate(components))
            #print len(resolved_comps_ids), len(components), len(filtered)
            #print filter(lambda (i, _): i in resolved_comps_ids, enumerate(components))
            components = map(lambda p: p[1], filtered)
            #print components

        for c in result:
            assert self.in_assembly(c.start)
            assert self.in_assembly(c.end)
            connections[-c.start] = Connection(-c.start, c.end, c.distance)
            connections[-c.end] = Connection(-c.end, c.start, c.distance)

        print "connections infered:", len(connections) / 2
        print "still {0} unresolved components".format(len(components))
        #print components
        return connections, components

