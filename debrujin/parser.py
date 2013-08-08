from utils import rc

#parse graph vertices


def graph_vertices(filename):
    graph = open(filename)
    ignorelines = 3
    for _ in range(ignorelines):
        graph.readline()

    while True:
        line = graph.readline()
        graph.readline()
        if '->' in line:
            break
        vinf = line.strip().split()
        v = vinf[0].strip('"')
        vid = int(v[:-1])
        length = int(vinf[1][3:])
        cov = int(vinf[2][2:-1])
        vid_conj = -1 * vid
        yield vid, vid_conj, length, cov

    graph.close()

#parse graph edges


def graph_edges(filename):
    graph = open(filename)
    for _ in range(3):
        line = graph.readline()
    defaultoverlap = int(line.strip().split()[1][3:-1])

    while '->' not in line:
        line = graph.readline()
    while True:
        edgeinfo = line.strip().split()
        v1id = int(edgeinfo[0][1:-2])
        if edgeinfo[0][-2] == '-':
            v1id = -1 * v1id
        v2id = int(edgeinfo[2][1:-2])
        if edgeinfo[2][-2] == '-':
            v2id = -1 * v2id

        newoverlap = defaultoverlap
        if len(edgeinfo) == 4:
            newoverlap = int(edgeinfo[3][3:-1])
        yield v1id, v2id, newoverlap
        line = graph.readline()
        if '}' in line:
            break
    graph.close()


def contigs_sequence(filename):
    contigs = open(filename)
    while True:
        idinf = contigs.readline().strip()
        if len(idinf) < 1:
            break
        seq = contigs.readline().strip()
        contigid = int(idinf.strip().split()[0][1:])
        yield contigid, seq
        yield -1 * contigid, rc(seq)



#for i in  graph_edges("data/spalgae-contigs.dot"):
#    print i
#for i in graph_vertices("data/spalgae-contigs.dot"):
#    print i
