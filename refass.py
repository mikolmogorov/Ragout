#!/usr/bin/env python

import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


import scripts.infer_connections as ic
import scripts.breakpoint_graph as bg
import scripts.sibelia_parser as sp


class DummyContig(sp.Contig):
    def __init__(self, length):
        sp.Contig.__init__(self, "")
        self.length = length

class Scaffold:
    def __init__(self, left, right, contigs):
        self.left = left
        self.right = right
        self.contigs = contigs


def calc_distance(offset, block_distance):
    if block_distance:
        return max(0, block_distance - offset)
    else:
        return 0


def extend_scaffolds(connections, sibelia_output):
    contigs = sibelia_output.contigs
    contig_index = sibelia_output.build_contig_index()
    scaffolds = []
    visited = set()

    def extend_scaffold(contig):
        visited.add(contig)
        scf = Scaffold(contig.blocks[0], contig.blocks[-1], [contig])
        scaffolds.append(scf)

        #go right
        while scf.right in connections:
            adjacent = connections[scf.right].end
            block_distance = connections[scf.right].distance
            assert len(contig_index[abs(adjacent)]) == 1

            contig = contigs[contig_index[abs(adjacent)][0]]
            if contig in visited:
                break

            if contig.blocks[0] == adjacent:
                offset = sibelia_output.block_offset(abs(adjacent), contig.name)
                reverse = scf.contigs[-1].sign > 0
                offset += sibelia_output.block_offset(abs(scf.right), scf.contigs[-1].name, reverse)
                #distance = max(0, block_distance - offset)
                distance = calc_distance(offset, block_distance)

                scf.contigs.append(DummyContig(distance))
                scf.contigs.append(contig)
                scf.right = contig.blocks[-1]
                visited.add(contig)
                continue

            if -contig.blocks[-1] == adjacent:
                offset = sibelia_output.block_offset(abs(adjacent), contig.name, True)
                reverse = scf.contigs[-1].sign > 0
                offset += sibelia_output.block_offset(abs(scf.right), scf.contigs[-1].name, reverse)
                distance = calc_distance(offset, block_distance)

                scf.contigs.append(DummyContig(distance))
                scf.contigs.append(contig)
                scf.contigs[-1].sign = -1
                scf.right = -contig.blocks[0]
                visited.add(contig)
                continue

            break

        #go left
        while -scf.left in connections:
            adjacent = -connections[-scf.left].end
            block_distance = connections[-scf.left].distance
            assert len(contig_index[abs(adjacent)]) == 1

            contig = contigs[contig_index[abs(adjacent)][0]]
            if contig in visited:
                break

            if contig.blocks[-1] == adjacent:
                offset = sibelia_output.block_offset(abs(adjacent), contig.name, True)
                reverse = scf.contigs[0].sign < 0
                offset += sibelia_output.block_offset(abs(scf.left), scf.contigs[0].name, reverse)
                #distance = max(0, block_distance - offset)
                distance = calc_distance(offset, block_distance)

                scf.contigs.insert(0, DummyContig(distance))
                scf.contigs.insert(0, contig)
                scf.left = contig.blocks[0]
                visited.add(contig)
                continue

            if -contig.blocks[0] == adjacent:
                offset = sibelia_output.block_offset(abs(adjacent), contig.name)
                reverse = scf.contigs[0].sign < 0
                offset += sibelia_output.block_offset(abs(scf.left), scf.contigs[0].name, reverse)
                #distance = max(0, block_distance - offset)
                distance = calc_distance(offset, block_distance)

                scf.contigs.insert(0, DummyContig(distance))
                scf.contigs.insert(0, contig)
                scf.contigs[0].sign = -1
                scf.left = -contig.blocks[-1]
                visited.add(contig)
                continue

            break

    for contig in contigs:
        if contig not in visited:
            extend_scaffold(contig)
    return scaffolds


def get_scaffolds(connections, sibelia_output):
    contigs = sibelia_output.contigs
    contig_index = sibelia_output.build_contig_index()
    scaffolds = extend_scaffolds(connections, sibelia_output)
    scaffolds = filter(lambda s: len(s.contigs) > 1, scaffolds)
    for scf in scaffolds:
        contigs = filter(lambda c : not isinstance(c, DummyContig), scf.contigs)
        for contig in contigs:
            if contig.sign > 0:
                print contig.blocks,
            else:
                print map(lambda b: -b, contig.blocks)[::-1],
        print ""
    return scaffolds


def parse_contigs(contigs_file):
    contigs = SeqIO.parse(contigs_file, "fasta")
    seqs = [seq for seq in contigs]
    names = [contig.id for contig in seqs]
    return seqs, names


def output_scaffolds(contigs_seqs, scaffolds, out_scaffolds, out_order, write_contigs=False):
    MIN_CONTIG_LEN = 500
    out_stream = open(out_scaffolds, "w")
    if out_order:
        out_order_stream = open(out_order, "w")
    queue = {}
    for rec in contigs_seqs:
        queue[rec.id] = rec.seq

    counter = 0
    for scf in scaffolds:
        name = "scaffold{0}".format(counter)
        counter += 1
        scf_seq = Seq("")
        buffer = ""
        if out_order:
            out_order_stream.write(">" + name + "\n")

        for i, contig in enumerate(scf.contigs):
            if isinstance(contig, DummyContig):
                buffer = "N" * max(11, contig.length) #for quasts`s scaffold splitting
                continue

            cont_seq = queue[contig.name]
            del queue[contig.name]

            if contig.sign < 0:
                cont_seq = cont_seq.reverse_complement()

            if i > 0:
                #check for overlapping
                overlap = False
                for window in xrange(5, 100):
                    if str(scf_seq)[-window:] == str(cont_seq)[0:window]:
                        #assert overlap == False
                        cont_seq = cont_seq[window:]
                        overlap = True
                if not overlap:
                    scf_seq += Seq(buffer)
                    if out_order:
                        out_order_stream.write("gaps {0}\n".format(len(buffer)))
            buffer = ""
            scf_seq += cont_seq
            if out_order:
                sign = "+" if contig.sign > 0 else "-"
                out_order_stream.write(sign + contig.name + "\n")

        SeqIO.write(SeqRecord(scf_seq, id=name, description=""), out_stream, "fasta")

    count = 0
    for h, seq in queue.iteritems():
        if len(seq) > MIN_CONTIG_LEN:
            count += 1
    print "Done,", count, "contigs left"

    counter = 0
    if write_contigs:
        for i, seq in enumerate(queue.itervalues()):
            SeqIO.write(SeqRecord(seq, id="contig{0}".format(i), description=""),
                                out_stream, "fasta")


def do_job(sibelia_dir, contigs_file, out_scaffolds, out_order, out_graph):
    contigs_seqs, contig_names = parse_contigs(contigs_file)
    sibelia_output = sp.SibeliaOutput(sibelia_dir, contig_names)

    graph = bg.build_graph(sibelia_output)

    adj_finder = ic.AdjacencyFinder(graph, sibelia_output)
    connections, unresolved_components = adj_finder.find_adjacencies()
    scaffolds = get_scaffolds(connections, sibelia_output)

    output_scaffolds(contigs_seqs, scaffolds, out_scaffolds, out_order, False)
    if out_graph:
        bg.output_graph(graph, unresolved_components, open(out_graph, "w"))


def main():
    parser = argparse.ArgumentParser(description="Tool for reference-assisted assembly")
    parser.add_argument("-s", action="store", metavar="sibelia_dir", dest="sibelia_dir",
                        required=True, help="Directory with Sibelia output")
    parser.add_argument("-c", action="store", metavar="contigs_file", dest="contigs_file",
                        required=True, help="File with contigs in fasta format")
    parser.add_argument("-o", action="store", metavar="output_file", dest="output_file",
                        default="scaffolds.fasta", help="Output scaffolds file (default: scaffolds.fasta)")
    parser.add_argument("-g", action="store", metavar="graph_file", dest="graph_file",
                        default=None, help="Output file for breakpoint graph (default: Not set)")
    parser.add_argument("-r", action="store", metavar="order_file", dest="order_file",
                        default=None, help="Output file for contigs order (default: Not set)")

    args = parser.parse_args()
    do_job(args.sibelia_dir, args.contigs_file, args.output_file, args.order_file, args.graph_file)

if __name__ == "__main__":
    main()
