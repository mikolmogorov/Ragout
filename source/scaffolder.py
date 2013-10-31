from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import sibelia_parser as sp
from datatypes import *


def extend_scaffolds(connections, sibelia_output):
    contigs = sibelia_output.get_filtered_contigs()
    contig_index = sp.build_contig_index(contigs)

    scaffolds = []
    visited = set()
    counter = [0]

    def extend_scaffold(contig):
        visited.add(contig)
        scf_name = "scaffold{0}".format(counter[0])
        counter[0] += 1
        scf = Scaffold.with_contigs(scf_name, contig.blocks[0], contig.blocks[-1], [contig])
        scaffolds.append(scf)

        #go right
        while scf.right in connections:
            adjacent = connections[scf.right].end

            assert len(contig_index[abs(adjacent)]) == 1

            contig = contigs[contig_index[abs(adjacent)][0]]
            if contig in visited:
                break

            if contig.blocks[0] == adjacent:
                scf.contigs.append(contig)
                scf.right = contig.blocks[-1]
                visited.add(contig)
                continue

            if -contig.blocks[-1] == adjacent:
                scf.contigs.append(contig)
                scf.contigs[-1].sign = -1
                scf.right = -contig.blocks[0]
                visited.add(contig)
                continue

            break

        #go left
        while -scf.left in connections:
            adjacent = -connections[-scf.left].end

            assert len(contig_index[abs(adjacent)]) == 1

            contig = contigs[contig_index[abs(adjacent)][0]]
            if contig in visited:
                break

            if contig.blocks[-1] == adjacent:
                scf.contigs.insert(0, contig)
                scf.left = contig.blocks[0]
                visited.add(contig)
                continue

            if -contig.blocks[0] == adjacent:
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
    scaffolds = extend_scaffolds(connections, sibelia_output)
    scaffolds = filter(lambda s: len(s.contigs) > 1, scaffolds)
    print "Done, {0} scaffolds".format(len(scaffolds))
    return scaffolds


def output_scaffolds(contigs_fasta, scaffolds, out_fasta, out_order):
    MIN_CONTIG_LEN = 0

    out_fasta_stream = open(out_fasta, "w")
    if out_order:
        out_order_stream = open(out_order, "w")
    used_contigs = set()

    for scf in scaffolds:
        scf_seq = Seq("")
        buffer = ""
        out_order_stream.write(">" + scf.name + "\n")

        for i, contig in enumerate(scf.contigs):
            cont_seq = contigs_fasta[contig.name]
            used_contigs.add(contig.name)

            if contig.sign < 0:
                cont_seq = cont_seq.reverse_complement()

            scf_seq += Seq("N" * 11)
            scf_seq += cont_seq

            out_order_stream.write(str(contig) + "\n")

        SeqIO.write(SeqRecord(scf_seq, id=scf.name, description=""), out_fasta_stream, "fasta")

    count = 0
    for h, seq in contigs_fasta.iteritems():
        if len(seq) > MIN_CONTIG_LEN and h not in used_contigs:
            count += 1
    print "Done,", count, "of", len(contigs_fasta), "contigs left,"
