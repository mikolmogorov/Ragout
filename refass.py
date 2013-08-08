#!/usr/bin/env python

import sys
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


import scripts.infer_connections as ic
import scripts.breakpoint_graph as bg
import scripts.sibelia_parser as sp


class DumbContig(sp.Contig):
    def __init__(self, length):
        sp.Contig.__init__(self, "")
        self.length = length

class Scaffold:
    def __init__(self, left, right, contigs):
        self.left = left
        self.right = right
        self.contigs = contigs

#def sign_of(val):
#    return math.copysign(1, val)


def extend_scaffolds(contigs, contig_index, connections):
    scaffolds = []
    visited = set()

    def extend_scaffold(contig):
        visited.add(contig)
        scf = Scaffold(contig.blocks[0], contig.blocks[-1], [contig])
        scaffolds.append(scf)

        #go right
        while scf.right in connections:
            adjacent = connections[scf.right].end
            distance = connections[scf.right].distance
            assert len(contig_index[abs(adjacent)]) == 1

            contig = contigs[contig_index[abs(adjacent)][0]]
            if contig in visited:
                break

            if contig.blocks[0] == adjacent:
                scf.contigs.append(DumbContig(distance))
                scf.contigs.append(contig)
                scf.right = contig.blocks[-1]
                visited.add(contig)
                continue

            if -contig.blocks[-1] == adjacent:
                scf.contigs.append(DumbContig(distance))
                scf.contigs.append(contig)
                scf.contigs[-1].sign = -1
                scf.right = -contig.blocks[0]
                visited.add(contig)
                continue

            break

        #go left
        while -scf.left in connections:
            adjacent = -connections[-scf.left].end
            distance = connections[-scf.left].distance
            assert len(contig_index[abs(adjacent)]) == 1

            contig = contigs[contig_index[abs(adjacent)][0]]
            if contig in visited:
                break

            if contig.blocks[-1] == adjacent:
                scf.contigs.insert(0, DumbContig(distance))
                scf.contigs.insert(0, contig)
                scf.left = contig.blocks[0]
                visited.add(contig)
                continue

            if -contig.blocks[0] == adjacent:
                scf.contigs.insert(0, DumbContig(distance))
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


def get_scaffolds(contigs, contig_index, connections):
    scaffolds = extend_scaffolds(contigs, contig_index, connections)
    scaffolds = filter(lambda s: len(s.contigs) > 1, scaffolds)
    for scf in scaffolds:
        contigs = filter(lambda c : not isinstance(c, DumbContig), scf.contigs)
        for contig in contigs:
            if contig.sign > 0:
                print contig.blocks,
            else:
                print map(lambda b: -b, contig.blocks)[::-1],
        print ""
    return scaffolds


def output_scaffolds(input_contigs, scaffolds, out_file, write_contigs=False):
    MIN_CONTIG_LEN = 500
    contigs = SeqIO.parse(input_contigs, "fasta")
    out_stream = open(out_file, "w")
    queue = {}
    for rec in contigs:
        queue[rec.id] = rec.seq

    counter = 0
    for scf in scaffolds:
        scf_seq = Seq("")
        buffer = ""

        for i, contig in enumerate(scf.contigs):
            if isinstance(contig, DumbContig):
                buffer = "N" * contig.length
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
                        assert overlap == False
                        cont_seq = cont_seq[window:]
                        overlap = True
                if not overlap:
                    scf_seq += Seq(buffer)
                #scf_seq += Seq("N" * 11)    #easy mode
            buffer = ""
            scf_seq += cont_seq

        name = "scaffold{0}".format(counter)
        counter += 1
        SeqIO.write(SeqRecord(scf_seq, id=name, description=""), out_stream, "fasta")

    count = 0
    for h, seq in queue.iteritems():
        if len(seq) > MIN_CONTIG_LEN:
            count += 1
    print "Done,", count, "contigs left"
    if write_contigs:
        for h, seq in queue.iteritems():
            SeqIO.write(SeqRecord(seq, id=h, description=""), out_stream, "fasta")


def do_job(sibelia_dir, contigs_file, out_scaffolds, out_graph):
    coords_file = os.path.join(sibelia_dir, "blocks_coords.txt")
    permutations_file = os.path.join(sibelia_dir, "genomes_permutations.txt")

    blocks_coords, seqid = sp.parse_coords_file(coords_file)
    permutations, contigs = sp.parse_permutations_file(permutations_file)
    contig_index = sp.build_contig_index(contigs)
    num_references = len(permutations)

    graph = bg.build_graph(permutations, blocks_coords)
    conected_comps = bg.get_connected_components(graph)
    connections = ic.simple_connections(graph, conected_comps, contigs, contig_index, num_references)
    scaffolds = get_scaffolds(contigs, contig_index, connections)

    output_scaffolds(contigs_file, scaffolds, out_scaffolds)
    if out_graph:
        bg.output_graph(graph, open(out_graph, "w"))


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
    args = parser.parse_args()
    do_job(args.sibelia_dir, args.contigs_file, args.output_file, args.graph_file)

if __name__ == "__main__":
    main()
