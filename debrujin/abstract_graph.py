import logging


class Abstract_Vertex(object):
    def __init__(self, vid):
        self.inne = []  # incoming edges
        self.oute = []  # outgoing edges
        self.vid = vid  # vertexid

    def outv(self):
        return [e.v2 for e in self.oute]

    def innv(self):
        return [e.v1 for e in self.inne]

    def __hash__(self):
        return self.vid


class Abstract_Edge(object):
    def __init__(self, v1, v2, eid):
        self.v1, self.v2 = v1, v2
        v1.oute.append(self)
        v2.inne.append(self)
        self.eid = eid

    def __hash__(self):
        return self.eid

    def length(self):
        pass


class Abstract_Graph(object):

    def __init__(self):
        self.es = {}  # key -->edges
        self.vs = {}  # key -->vertices
        self.logger = logging.getLogger('Refass')
