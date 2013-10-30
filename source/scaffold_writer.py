from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


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
