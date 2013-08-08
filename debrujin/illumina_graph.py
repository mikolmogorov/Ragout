from abstract_graph import Abstract_Graph, Abstract_Edge, Abstract_Vertex
import logging
import parser
import sys
from utils import conjugate
from sets import Set
import math



class OVertex(Abstract_Vertex):
    def __init__(self, vid, conj, vlen, coverage):
        Abstract_Vertex.__init__(self, vid)
        conjugate(self, conj)
        self.length = vlen
        self.seq = None
        self.cov = coverage
        self.positions_on_ref = []

    def __repr__(self):
        return '"%s%s" [l=%d C=%d]' % (self.Name(), self.Strand(),
                                       self.length, self.cov)

    def Strand(self):
        if self.vid > 0:
            return '+'
        if self.vid < 0:
            return '-'

    def Name(self):
        return str(abs(self.vid))


class OEdge(Abstract_Edge):
    def __init__(self, edgeid, v1, v2, overlap):
        Abstract_Edge.__init__(self, v1, v2, edgeid)
        self.ovl = overlap
        self.conj = None

    def __repr__(self):
        return '"%s%s" -> "%s%s" [d=%d]' % \
            (self.v1.Name(), self.v1.Strand(),
             self.v2.Name(), self.v2.Strand(), self.ovl)


class Graph(Abstract_Graph):
    def __init__(self):
        Abstract_Graph.__init__(self)
        self.max_eid = 0
        self.logger = logging.getLogger('Algae')

    def add_vertex(self, vid, conj_id, vlen, cov):
        if vid in self.vs:
            v = self.vs[vid]
            if v.conj.vid != conj_id:
                sys.exit("Inconsistent conjugate vertex")
            if v.length != vlen:
                sys.exit("Inconsistent vertex length")
            if v.cov != cov:
                sys.exit("Inconsistent vertex coverage")
            return v
        conj = self.vs.get(conj_id, None)
        v = OVertex(vid, conj, vlen, cov)
        self.vs[vid] = v
        return v

    def add_edge(self, v1id, v2id, overlap, eid=None):
        v1 = self.vs[v1id]
        v2 = self.vs[v2id]
        if eid is not None and eid in self.es:
            if (self.es[eid].v1.vid != v1.vid and
                self.es[eid].v2.vid != v2.vid):
                sys.exit('Inconsistent eid')
            else:
                return self.es[eid]

        if v2 in v1.outv():
            for e in v1.oute:
                if e.v2 == v2:
                    if (eid is None or e.eid == eid):
                        return e
                    sys.exit('Inconsistent eid')
            sys.exit('Edge not found!')

        if eid is None:
            while self.max_eid in self.es.keys():
                self.max_eid += 1
            neweid = self.max_eid
            self.max_eid += 1
        else:
            if eid in self.es.keys():
                sys.exit("EID exists")
            neweid = eid
        e = OEdge(neweid, v1, v2, overlap)
        v1 = self.vs[v1id]
        v2 = self.vs[v2id]
        self.es[e.eid] = e
        #if the conjugate already exists, add conjugate
        v1conj = v1.conj
        v2conj = v2.conj
        if v1conj in v2conj.outv():
            ej = None
            for ec in v2conj.oute:
                if ec.v2 == v1conj:
                    ej = ec
            if not ej:
                sys.exit('Conjugate edge not found!')
            conjugate(e, ej)
        return e

    def remove_edge(self, e):
#         if e not in self.vs[e.v1.vid].oute:
#             sys.exit("Edge not in oute")
#         if e not in self.vs[e.v2.vid].oute:
#             sys.exit("Edge not in inne")
        for e1 in self.vs[e.v1.vid].oute:
            if e1.eid == e.eid:
                self.vs[e.v1.vid].oute.remove(e1)
                break
        for e2 in self.vs[e.v2.vid].inne:
            if e2.eid == e.eid:
                self.vs[e.v2.vid].inne.remove(e2)
                break
        self.es.pop(e.eid)
        if e.conj is not None:
            for e1c in self.vs[e.conj.v1.vid].oute:
                if e1c.eid == e.conj.eid:
                    self.vs[e.conj.v1.vid].oute.remove(e1c)
                    break
            for e2c in self.vs[e.conj.v2.vid].inne:
                if e2c.eid == e.conj.eid:
                    self.vs[e.conj.v2.vid].inne.remove(e2c)
                    break
            self.es.pop(e.conj.eid)
        return None

    def has_edge(self, v1id, v2id, overlap, offset_margin=0):
        if v1id not in self.vs or v2id not in self.vs:
            return False
        if v2id not in self.vs[v1id].outv():
            return False
        for e in self.vs[v1id].oute:
            if e.v2.vid == v2id and (overlap >= (1 - offset_margin) * e.ovl and
                                     overlap <= (1 + offset_margin) * e.ovl):
                return True
        return False

    def find_path(self, v1id, v2id, min_offset, max_offset, depth=10):
        if max_offset < 0:
            return None
        if depth < 0:
            return None
        v1 = self.vs[v1id]
        if v1id == v2id:
            if min_offset <= 0 and max_offset >= 0:
                return []
        for ve in v1.oute:
            off = v1.length + ve.ovl
            p = self.find_path(ve.v2.vid, v2id, min_offset - off,
                               max_offset - off, depth - 1)
            if p is not None:
                return [ve] + p
        return None

    def add_seq(self, vid, seq):
        self.vs[vid].seq = seq

    def load_dna(self, filename):
        for contigid, seq in parser.contigs_sequence(filename):
            self.add_seq(contigid, seq)

    def add_cov(self, vid, cvr):
        self.vs[vid].cvr = cvr

    def path_seq(self, path, offset_margin=0.01, fill_gaps=True):
        vertices = []
        overlapping = []
        for e in path:
            vertices.append(self.vs[e.v1.vid])
            overlapping.append(e.ovl)
        if len(path) == 0:
            return ''
        vertices.append(self.vs[path[-1].v2.vid])

        seq = vertices[0].seq
        for i in range(1, len(vertices)):
            e = path[i - 1]
            v1 = self.vs[e.v1.vid]
            v2 = self.vs[e.v2.vid]
            if path[i - 1].ovl > 0:
                if fill_gaps == False:
                    seq += 'N' * path[i - 1].ovl
                    seq += vertices[i].seq
                else:
                    offset = v1.length + e.ovl
                    min_offset = (1 - offset_margin) * offset
                    max_offset = (1 + offset_margin) * offset
    #                print v1.vid, e.ovl, v2.vid, min_offset, max_offset
                    subpath = self.find_path(v1.vid, v2.vid,
                                             min_offset, max_offset)
                    if subpath is not None:
                        if len(subpath) == 1:
                            seq += 'N' * path[i - 1].ovl
                            seq += vertices[i].seq
                        else:
                            subpathseq = self.path_seq(subpath, offset_margin)
                            seq += subpathseq[v1.length:]
                    else:
                        seq += 'N' * path[i - 1].ovl
                        seq += vertices[i].seq
            else:
                if vertices[i].seq is None:
                    print i
                seq += vertices[i].seq[-1 * overlapping[i - 1]:]
        return seq

    def load(self, graph_filename, seq_filename=None):
        for vid, conj, length, cov in parser.graph_vertices(graph_filename):
            self.add_vertex(vid, conj, length, cov)
            self.add_vertex(conj, vid, length, cov)

        for v1, v2, overlap in parser.graph_edges(graph_filename):
            self.add_edge(v1, v2, overlap)

        if seq_filename is not None:
            for vid, vseq in parser.contigs_sequence(seq_filename):
                if vid in self.vs.keys():
                    if self.vs[vid].length != len(vseq):
                        print('Graph.load: Vertex length not consistent' +
                                 '\t'.join(map(str, [vid, self.vs[vid].length,
                                               len(vseq)])))
                        if len(vseq) < self.vs[vid].length:
                            if vid > 0:
                                self.add_seq(vid, vseq +
                                             ('N' *
                                              (self.vs[vid].length -
                                               len(vseq))))
                            else:
                                self.add_seq(vid, ('N' * (self.vs[vid].length -
                                               len(vseq))) + vseq)
                        else:
                            if vid > 0:
                                self.add_seq(vid, vseq[:self.vs[vid].length])
                            else:
                                self.add_seq(vid,
                                             vseq[-1 * self.vs[vid].length:])
                        continue
                    self.add_seq(vid, vseq)
                else:
                    length = len(vseq)
                    self.add_vertex(vid, -1 * vid, length, 0)
                    self.add_seq(vid, vseq)

    def get_vid(self, name, strand):
        if strand == '+':
            return int(name)
        else:
            return -1 * int(name)
    """
    vid: vertex id
    positions:  if positions < 0 it means in the reverse complement

    """
    def add_positions_ref(self, vid, positions):
        self.vs[vid].positions_on_ref = positions[:]

    def __repr__(self):
        S = "digraph adj {\n" + "graph [k=64]\n" + "edge [d=-63]\n"
        vsvalues = self.vs.values()[:]
        vsvalues.sort(lambda x, y: abs(x.vid) - abs(y.vid))
        for v in vsvalues:
            S = S + v.__repr__() + '\n'
        esvalues = self.es.values()[:]
        esvalues.sort(lambda x, y: abs(x.v1.vid) - abs(y.v1.vid))
        for e in esvalues:
            S = S + e.__repr__() + '\n'
        S = S + '}'
        return S

    def remove_transitive_edges(self, offset_margin):
        edges_removed = []
        edgeIDs = self.es.keys()[:]
        for eid in edgeIDs:
            if eid not in self.es.keys():
                continue
            e = self.es[eid]
#             if e.ovl < 0:
#                 continue
            for enext in self.vs[e.v1.vid].oute:
                if e.eid == enext.eid:
                    continue
                vnext = enext.v2
                offset = e.v1.length + e.ovl
                min_offset = ((1 - offset_margin) * offset -
                              (enext.v1.length + enext.ovl))
                max_offset = ((1 + offset_margin) * offset -
                              (enext.v1.length + enext.ovl))
                p = self.find_path(vnext.vid, e.v2.vid, min_offset, max_offset)
                if p is not None:
                    self.remove_edge(e)
                    edges_removed.append(e)
                    if e.conj is not None:
                        edges_removed.append(e.conj)
                    break
        return edges_removed

    def get_simple_paths(self):
        vid_seen = Set([])
        paths = []
        for vid in self.vs:
            if vid in vid_seen:
                continue
            v = self.vs[vid]
            if len(v.inne) > 1 or len(v.oute) > 1:
                continue
            if len(v.oute) <= 1:
                if len(v.oute) == 0 and len(v.inne) == 0:
                    continue
                if len(v.inne) == 1 and (len(v.inne[0].v1.oute) <= 1 and
                                         len(v.inne[0].v1.inne) <= 1):
                    continue
                vid_seen.add(v.vid)
                vid_seen.add(v.conj.vid)
                if len(v.inne) == 1:
                    p = [v.inne[0]]
                else:
                    p = []
                if len(v.oute) == 1:
                    ce = v.oute[0]
                    p.append(ce)
                    vid_seen.add(ce.v2.vid)
                    vid_seen.add(ce.v2.conj.vid)
                    while len(ce.v2.oute) == 1 and len(ce.v2.inne) == 1:
                        ce = ce.v2.oute[0]
                        p.append(ce)
                        vid_seen.add(ce.v2.vid)
                        vid_seen.add(ce.v2.conj.vid)

                paths.append(p)
        return paths

    def Subgraph(self, vID, radius):
        queue = []
        visited = []

        #Do BFS recurslly with input radius wished
        def BFS(current, radius):
            if radius < 0:
                return
            else:
                if len(current.oute) > 0:
                    queue.append(current.oute)
                    visited.append(current.vid)
                for path in current.oute:
                    #check for repeats
                    if path.v2.vid not in visited:
                        BFS(path.v2, radius - 1)

        BFS(self.vs[vID], radius - 1)
        return queue

    def get_seq(self, contigs_list):

        if len(contigs_list) < 1:
            return "Non valid"

        seq = ''
        ovlap = 0

        #add the first seq
        seq += self.vs[contigs_list[0]].seq
        #add the rest and - the ovlap
        for i in xrange(1, len(contigs_list)):
            # get overlap value
            for x in self.vs[contigs_list[i - 1]].oute:
                if abs(contigs_list[i]) == abs(x.v2.vid):
                    ovlap = -1 * x.ovl
            seq += self.vs[contigs_list[i]].seq[(ovlap + 1)::]

        return "".join(seq.split())

#g = Graph()
#g.load("data/spalgae-contigs.dot",None)
#branches = 0
#small_branches = 0
#
#for i in g.vs.values():
#    if len(i.outv()) > 1 and len(i.innv())> 1:
#
##        print "multiple"
#        print i.length
#        if i.length < 100:
#            small_branches +=1
#        branches +=1
##        print "Middle", i
##        print "Out" , i.outv()
##        print "Inn", i.innv()
#portion = "small contigs", small_branches*1.0/branches
#print "small contigs portion", portion

#g.load_dna("data/spalgae-contigs.fasta")
#print g.path_seq([g.es[0]])

#for i in g.vs.values():
#    print i

#for i in g.es.values():
#    print i
