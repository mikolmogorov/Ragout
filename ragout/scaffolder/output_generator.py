#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

from itertools import repeat
import logging

from ragout.parsers.fasta_parser import write_fasta_dict, reverse_complement

logger = logging.getLogger()

MIN_GAP = 11

def output_links(scaffolds, out_links):
    """
    Outputs pretty table with information about adjacencies
    """
    HEADER = ["contig_1", "contig_2", "gap", "ref_support", "~>"]
    COL_GAP = 4

    with open(out_links, "w") as f:
        for scf in scaffolds:
            rows = []
            for left, right in zip(scf.contigs[:-1], scf.contigs[1:]):
                supp_genomes = ",".join(sorted(left.link.supporting_genomes))
                supp_assembly = "*" if left.link.supporting_assembly else " "
                rows.append([str(left), str(right), str(left.link.gap),
                            supp_genomes, supp_assembly])

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
    total_contigs = 0
    total_len = 0
    #overlap = 0
    for scf in scaffolds:
        scf_seqs = []
        for contig in scf.contigs:
            ###
            sep = contig.seq_name.find("[")
            if sep == -1:
                seq_name = contig.seq_name
                cont_seq = contigs_fasta[contig.seq_name]
            else:
                seq_name = contig.seq_name[:sep]
                seg_start, seg_end = \
                        map(int, contig.seq_name[sep + 1 : -1].split(":"))
                cont_seq = contigs_fasta[seq_name][seg_start:seg_end]
            ###

            used_contigs.add(seq_name)
            total_contigs += 1
            total_len += len(cont_seq)

            if contig.sign < 0:
                cont_seq = reverse_complement(cont_seq)

            if contig.link.gap >= 0:
                scf_seqs.append(cont_seq)
                scf_seqs.append("N" * (max(MIN_GAP, contig.link.gap)))
            else:
                scf_seqs.append(cont_seq[:contig.link.gap])


        scf_seq = "".join(scf_seqs)
        scf_length.append(len(scf_seq))
        out_fasta_dict[scf.name] = scf_seq
    write_fasta_dict(out_fasta_dict, out_file)

    #add some statistics
    used_unique = 0
    used_len = 0
    unused_count = 0
    unused_len = 0
    for h in contigs_fasta:
        if h in used_contigs:
            used_unique += 1
            used_len += len(contigs_fasta[h])
        else:
            unused_count += 1
            unused_len += len(contigs_fasta[h])
    assembly_len = unused_len + used_len
    used_perc = 100 * float(used_len) / assembly_len
    unused_perc = 100 * float(unused_len) / assembly_len
    contigs_length = [len(c) for c in contigs_fasta.values()]

    logger.info("Assembly statistics:\n\n"
                "\tScaffolds:\t\t{0}\n"
                "\tUnique contigs:\t\t{1}\n"
                "\tUnique contigs length:\t{2} ({3:2.4}%)\n"
                "\tTotal contigs:\t\t{4}\n"
                "\tTotal contigs length:\t{5}\n"
                "\tUnused contigs count:\t{6}\n"
                "\tUnused contigs length:\t{7} ({8:2.4}%)\n"
                "\tContigs N50: \t\t{9}\n"
                "\tScaffolds N50:\t\t{10}\n"
                .format(len(scaffolds), used_unique, used_len, used_perc,
                        total_contigs, total_len, unused_count, unused_len,
                        unused_perc,
                        _calc_n50(contigs_length, unused_len + used_len),
                        _calc_n50(scf_length, unused_len + used_len)))


def _calc_n50(scaffolds_lengths, assembly_len):
    n50 = 0
    sum_len = 0
    for l in sorted(scaffolds_lengths, reverse=True):
        sum_len += l
        if sum_len > assembly_len / 2:
            n50 = l
            break
    return n50
