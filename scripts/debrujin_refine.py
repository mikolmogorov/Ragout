#!/usr/bin/env python

import networkx as nx
import sys
from collections import namedtuple
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


Scaffold = namedtuple("Scaffold", ["name", "contigs"])
Edge = namedtuple("Edge", ["start", "end"])


class Contig:
    def __init__(self, name, sign, gap = 0):
        self.name = name
        self.sign = sign
        self.gap = gap

    @staticmethod
    def from_sting(string):
        return Contig(*parse_contig_name(string))

    def __str__(self):
        sign = "+" if self.sign > 0 else "-"
        return sign + self.name


def load_graph(filename):
    graph = nx.DiGraph()
    edges = {}
    for line in open(filename, "r"):
        edge_id, start, end = line.strip().split(" ")

        graph.add_nodes_from([start, end])
        graph.add_edge(start, end)
        edges[edge_id] = Edge(start, end)
    return graph, edges


def parse_contig_name(string):
    if string[0] not in ["+", "-"]:
        return None

    sign = 1 if string[0] == "+" else -1
    return string[1:], sign


def parse_contigs_order(filename):
    scaffolds = []
    for line in open(filename, "r"):
        line = line.strip()
        if line.startswith(">"):
            scaffolds.append(Scaffold(line[1:], []))
        else:
            if line.startswith("gaps"):
                gaplen = int(line.split(" ")[1])
                scaffolds[-1].contigs[-1].gap = gaplen
            else:
                #name = line.strip("\n").replace("=", "_") #fix for quast
                #name, sign = parse_contig_(line.strip())
                scaffolds[-1].contigs.append(Contig.from_sting(line))
    return scaffolds


def insert_from_graph(graph_file, scaffolds_in):
    new_scaffolds = []
    graph, edges = load_graph(graph_file)

    for scf in scaffolds_in:
        new_scaffolds.append(Scaffold(scf.name, []))

        for prev_cont, new_cont in zip(scf.contigs[:-1], scf.contigs[1:]):
            src = edges[str(prev_cont)].end
            dst = edges[str(new_cont)].start

            new_scaffolds[-1].contigs.append(prev_cont)

            if src == dst:
                print "adjacent contigs {0} and {1}".format(prev_cont, new_cont)
                new_scaffolds[-1].contigs[-1].gaps = 0
                continue

            if not nx.has_path(graph, src, dst):
                print "no path between {0} and {1}".format(prev_cont, new_cont)
                continue

            paths = [p for p in nx.all_shortest_paths(graph, src, dst)]
            if len(paths) != 1:
                print "ambigious paths between {0} and {1}".format(prev_cont, new_cont)
                continue

            path = paths[0]
            if len(path) > 4:
                print "too long path between {0} and {1}".format(prev_cont, new_cont)
                continue

            #all is ok!
            print ("found path between {0} and {1} of length {2}"
                                    .format(prev_cont, new_cont, len(path)))
            for p_start, p_end in zip(paths[0][:-1], paths[0][1:]):
                #corresponging edge in graph
                found_edge = None
                for edge_id, edge in edges.iteritems():
                    if edge == Edge(p_start, p_end):
                        found_edge = edge_id
                        break
                assert found_edge
                #assert found_edge[1:] not in used_contigs

                new_scaffolds[-1].contigs[-1].gaps = 0
                new_scaffolds[-1].contigs.append(Contig.from_sting(found_edge))

        new_scaffolds[-1].contigs.append(new_cont)

    return new_scaffolds


def output_scaffolds(contigs_fasta, scaffolds, out_fasta, out_order, write_contigs=False):
    MIN_CONTIG_LEN = 500
    KMER = 55
    OVERLAP_DIST = [KMER]

    out_fasta_stream = open(out_fasta, "w")
    if out_order:
        out_order_stream = open(out_order, "w")
    used_contigs = set()

    gaps = 0
    for scf in scaffolds:
        scf_seq = Seq("")
        buffer = ""
        if out_order:
            out_order_stream.write(">" + scf.name + "\n")

        for i, contig in enumerate(scf.contigs):
            cont_seq = contigs_fasta[contig.name]

            if contig.sign < 0:
                cont_seq = cont_seq.reverse_complement()

            if gaps > 0:
                scf_seq += Seq("N" * max(11, gaps))
            elif i > 0:
                #check for overlapping
                overlap = False
                for window in OVERLAP_DIST:
                    if str(scf_seq)[-window:] == str(cont_seq)[0:window]:
                        cont_seq = cont_seq[window:]
                        overlap = True
                        #print "overlap!"
                if not overlap:
                    print "no overlap!"

            scf_seq += cont_seq
            used_contigs.add(contig.name)
            gaps = contig.gap

            if out_order:
                out_order_stream.write(str(contig) + "\ngaps {0}\n".format(contig.gap))

        SeqIO.write(SeqRecord(scf_seq, id=scf.name, description=""), out_fasta_stream, "fasta")

    count = 0
    for h, seq in contigs_fasta.iteritems():
        if len(seq) > MIN_CONTIG_LEN and h not in used_contigs:
            count += 1
    print "Done,", count, "contigs left"

    if write_contigs:
        for i, hdr in enumerate(contigs_fasta):
            if hdr in used_contigs:
                continue
            SeqIO.write(SeqRecord(contigs_fasta[hdr], id="contig{0}".format(i), description=""),
                                                                    out_fasta_stream, "fasta")


def refine_contigs(graph_file, scaffolds, contigs_in, scaffolds_out, order_out):
    graph, edges = load_graph(graph_file)
    contigs = {seq.id : seq.seq for seq in SeqIO.parse(contigs_in, "fasta")}

    new_scaffolds = insert_from_graph(graph_file, scaffolds)
    output_scaffolds(contigs, new_scaffolds, scaffolds_out, order_out, False)


def main():
    if len(sys.argv) < 3:
        print "Usage: debrujin spades_file contig_oder contig_file"
        return

    graph_filename, order_filename, contigs_filename = sys.argv[1], sys.argv[2], sys.argv[3]
    scaffolds = parse_contigs_order(order_filename)
    refine_contigs(graph_filename, scaffolds, contigs_filename, "new_scf.fasta", "new_order.txt")


if __name__ == "__main__":
    main()
