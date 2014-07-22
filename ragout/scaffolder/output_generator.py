#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

from itertools import repeat
import logging

from ragout.parsers.fasta_parser import write_fasta_dict, reverse_complement

logger = logging.getLogger()

def output_links(scaffolds, out_links):
    """
    Outputs pretty table with information about adjacencies
    """
    HEADER = ["contig_1", "contig_2", "gap", "~>", "ref_support"]
    COL_GAP = 4

    with open(out_links, "w") as f:
        for scf in scaffolds:
            rows = []
            for left, right in zip(scf.contigs[:-1], scf.contigs[1:]):
                supp_genomes = ",".join(sorted(left.link.supporting_genomes))
                supp_assembly = "*" if left.link.supporting_assembly else " "
                rows.append([str(left), str(right), str(left.link.gap),
                            supp_assembly, supp_genomes])

            col_widths = repeat(0)
            for row in [HEADER] + rows:
                col_widths = [max(len(v), w) for v, w in zip(row, col_widths)]
            line_len = sum(col_widths) + COL_GAP * len(col_widths)

            #header
            f.write("-" * line_len + "\n")
            f.write(scf.name + "\n")
            f.write("-" * line_len + "\n")
            for hdr, width in zip(HEADER, col_widths):
                f.write(hdr + (" " * (width - len(hdr) + COL_GAP)))
            f.write("\n" + "-" * line_len + "\n")

            #values
            for row in rows:
                for val, width in zip(row, col_widths):
                    f.write(val + (" " * (width - len(val) + COL_GAP)))
                f.write("\n")

            f.write("-" * line_len + "\n\n")


def output_fasta(contigs_fasta, scaffolds, out_file):
    """
    Outputs scaffodls to file in "fasta" format
    """
    logger.info("Generating FASTA output")
    used_contigs = set()
    out_fasta_dict = {}

    scf_length = []
    for scf in scaffolds:
        scf_seqs = []
        for contig in scf.contigs:
            cont_seq = contigs_fasta[contig.name]
            used_contigs.add(contig.name)

            if contig.sign < 0:
                cont_seq = reverse_complement(cont_seq)

            if contig.link.gap >= 0:
                scf_seqs.append(cont_seq)
                scf_seqs.append("N" * contig.link.gap)
            else:
                scf_seqs.append(cont_seq[:contig.link.gap])

        scf_seq = "".join(scf_seqs)
        scf_length.append(len(scf_seq))
        out_fasta_dict[scf.name] = scf_seq
    write_fasta_dict(out_fasta_dict, out_file)

    #add some statistics
    used_count = 0
    used_len = 0
    unused_count = 0
    unused_len = 0
    for h in contigs_fasta:
        if h in used_contigs:
            used_count += 1
            used_len += len(contigs_fasta[h])
        else:
            unused_count += 1
            unused_len += len(contigs_fasta[h])
    assembly_len = unused_len + used_len
    used_perc = 100 * float(used_len) / assembly_len
    unused_perc = 100 * float(unused_len) / assembly_len
    contigs_length = [len(c) for c in contigs_fasta.values()]

    logger.info("Assembly statistics:\n\n"
                "\tScaffolds count:\t{0}\n"
                "\tUsed contigs count:\t{1}\n"
                "\tUsed contigs length:\t{2} ({3:2.4}%)\n"
                "\tUnused contigs count:\t{4}\n"
                "\tUnused contigs length:\t{5} ({6:2.4}%)\n"
                "\tContigs N50: \t\t{7}\n"
                "\tScaffolds N50:\t\t{8}\n"
                .format(len(scaffolds), used_count, used_len, used_perc,
                        unused_count, unused_len, unused_perc,
                        _calc_n50(contigs_length, unused_len + used_len),
                        _calc_n50(scf_length, unused_len + used_len)))


#def output_order(scaffolds, out_order):
#    """
#    Outputs scaffolds to file in "ord" format
#    """
#    with open(out_order, "w") as f:
#        for scf in scaffolds:
#            f.write(">" + scf.name + "\n")
#            for contig in scf.contigs:
#                f.write(str(contig) + "\n")


def _calc_n50(scaffolds_lengths, assembly_len):
    n50 = 0
    sum_len = 0
    for l in sorted(scaffolds_lengths, reverse=True):
        sum_len += l
        if sum_len > assembly_len / 2:
            n50 = l
            break
    return n50
