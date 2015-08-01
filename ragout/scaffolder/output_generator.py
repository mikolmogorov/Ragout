#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

from itertools import repeat
import logging
import os

from ragout.parsers.fasta_parser import write_fasta_dict, reverse_complement
from ragout.__version__ import __version__

logger = logging.getLogger()

MIN_GAP = 11


def make_output(contigs, scaffolds, out_dir, out_prefix):
    out_links = os.path.join(out_dir, out_prefix + "_scaffolds.links")
    out_agp = os.path.join(out_dir, out_prefix + "_scaffolds.agp")
    out_chr = os.path.join(out_dir, out_prefix + "_scaffolds.fasta")
    out_unplaced = os.path.join(out_dir, out_prefix + "_unplaced.fasta")
    _fix_gaps(contigs, scaffolds)
    output_links(scaffolds, out_links)
    _output_agp(scaffolds, out_agp, out_prefix)
    _output_fasta(contigs, scaffolds, out_chr, out_unplaced)


def _fix_gaps(contigs, scaffolds):
    def get_seq(contig):
        seq_name, seg_start, seg_end = contig.name_with_coords()
        if seg_start is None:
            cont_seq = contigs[seq_name]
        else:
            cont_seq = contigs[seq_name][seg_start:seg_end]
        if contig.sign < 0:
            cont_seq = reverse_complement(cont_seq)
        return cont_seq

    def count_ns(cnt_1, cnt_2):
        seq_1, seq_2 = get_seq(cnt_1), get_seq(cnt_2)
        left_ns, right_ns = 0, 0
        for i in xrange(len(seq_1) - 1, 0, -1):
            if seq_1[i].upper() != "N":
                break
            left_ns += 1
        for i in xrange(len(seq_2) - 1):
            if seq_2[i].upper() != "N":
                break
            right_ns += 1
        return left_ns, right_ns

    for scf in scaffolds:
        for cnt_1, cnt_2 in zip(scf.contigs[:-1], scf.contigs[1:]):
            if cnt_1.link.supporting_assembly:
                cnt_1.trim_right(max(0, -cnt_1.link.gap))
                cnt_1.link.gap = max(0, cnt_1.link.gap)
                continue

            left_ns, right_ns = count_ns(cnt_1, cnt_2)

            cnt_1.trim_right(left_ns)
            cnt_2.trim_left(right_ns)
            cnt_1.link.gap += left_ns + right_ns
            cnt_1.link.gap = max(cnt_1.link.gap, MIN_GAP)


def _output_agp(scaffolds, out_agp, assembly_name):
    """
    Output file in NCBI AGP format
    """
    SHIFT = 1
    with open(out_agp, "w") as f:
        f.write("##agp-version  2.0\n")
        f.write("#ASSEMBLY NAME: {0}\n".format(assembly_name))
        f.write("#DESCRIPTION: Pseudochromosome assembly\n")
        f.write("#PROGRAM: Ragout v{0}\n".format(__version__))
        for scf in scaffolds:
            chr_pos = 0
            for contig_id, contig in enumerate(scf.contigs):
                chr_start = chr_pos
                chr_end = chr_pos + contig.length()
                chr_pos = chr_end + contig.link.gap
                cont_name, cont_start, cont_end = contig.name_with_coords()
                strand = "+" if contig.sign > 0 else "-"

                supp_genomes = sorted(contig.link.supporting_genomes)
                if contig.link.supporting_assembly:
                    supp_genomes.append("~>")
                support = ",".join(supp_genomes)

                contig_num = 2 * contig_id + 1
                gap_num = 2 * contig_id + 2
                cont_fields = [scf.name, chr_start + SHIFT, chr_end,
                               contig_num, "W", cont_name, cont_start + SHIFT,
                               cont_end, strand]
                f.write("\t".join(map(str, cont_fields)) + "\n")
                if contig.link.gap > 0:
                    gap_fields = [scf.name, chr_end + SHIFT, chr_pos, gap_num,
                                  "N", contig.link.gap, "scaffold", "yes", support]
                    f.write("\t".join(map(str, gap_fields)) + "\n")


def output_links(scaffolds, out_links):
    """
    Outputs pretty table with information about adjacencies
    """
    HEADER = ["sequence", "start", "length", "gap", "support"]
    COL_GAP = 4

    with open(out_links, "w") as f:
        for scf in scaffolds:
            rows = []
            cur_pos = 0

            for contig in scf.contigs:
                start = cur_pos
                cur_pos = start + contig.length() + contig.link.gap

                supp_genomes = sorted(contig.link.supporting_genomes)
                if contig.link.supporting_assembly:
                    supp_genomes.append("~>")
                support = ",".join(supp_genomes)

                rows.append([contig.signed_name(), str(start),
                            str(contig.length()), str(contig.link.gap),
                            support])

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


def _output_fasta(contigs_fasta, scaffolds, out_chr, out_unlocalized):
    """
    Outputs scaffodls to file in "fasta" format
    """
    logger.info("Generating FASTA output")
    used_contigs = set()
    out_chromosomes = {}

    scf_length = []
    total_contigs = 0
    total_len = 0
    for scf in scaffolds:
        scf_seqs = []
        #trim_left = 0
        for contig in scf.contigs:
            seq_name, seg_start, seg_end = contig.name_with_coords()
            if seg_start is None:
                cont_seq = contigs_fasta[seq_name]
            else:
                cont_seq = contigs_fasta[seq_name][seg_start:seg_end]
            if contig.sign < 0:
                cont_seq = reverse_complement(cont_seq)

            #if contig.link.trim_left > 0:
            #    scf_seqs.append(cont_seq[trim_left : -contig.link.trim_left])
            scf_seqs.append(cont_seq)
            if contig.link.gap > 0:
                scf_seqs.append("N" * contig.link.gap)
            #trim_left = contig.link.trim_right

            used_contigs.add(seq_name)
            total_len += len(cont_seq)

        total_contigs += len(scf.contigs)
        scf_seq = "".join(scf_seqs)
        scf_length.append(len(scf_seq))
        out_chromosomes[scf.name] = scf_seq
    write_fasta_dict(out_chromosomes, out_chr)

    #unused contigs
    unused_fasta = {}
    for cname in contigs_fasta:
        if cname not in used_contigs:
            unused_fasta[cname] = contigs_fasta[cname]
    write_fasta_dict(unused_fasta, out_unlocalized)

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
