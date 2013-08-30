from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def output_scaffolds(contigs_fasta, scaffolds, out_fasta, out_order, kmer_len, write_contigs=False):
    MIN_CONTIG_LEN = 0
    #KMER = kmer_len
    OVERLAP_DIST = [kmer_len]

    out_fasta_stream = open(out_fasta, "w")
    if out_order:
        out_order_stream = open(out_order, "w")
    used_contigs = set()

    gaps = 0
    nonoverlap = 0
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
                    scf_seq += Seq("N" * 11)
                    #print "no overlap!"
                    nonoverlap += 1

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
    print "Done,", count, "contigs left,", nonoverlap, "contigs were not overlapping"

    if write_contigs:
        for i, hdr in enumerate(contigs_fasta):
            if hdr in used_contigs:
                continue
            SeqIO.write(SeqRecord(contigs_fasta[hdr], id="contig{0}".format(i), description=""),
                                                                    out_fasta_stream, "fasta")

