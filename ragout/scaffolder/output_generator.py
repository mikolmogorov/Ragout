#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

from __future__ import absolute_import
from __future__ import division
from itertools import repeat
from collections import defaultdict
import logging
import os

from ragout.parsers.fasta_parser import write_fasta_dict, reverse_complement
from ragout.__version__ import __version__
import ragout.shared.config as config
from ragout.six.moves import map, range, zip

logger = logging.getLogger()


class OutputGenerator:
    def __init__(self, fragments_fasta, scaffolds):
        self.fragments_fasta = fragments_fasta
        self.scaffolds = scaffolds

        self.unplaced_fasta = None
        self.scaffolds_fasta = None
        self.used_fragments_len = None
        self.introduced_gap_len = None

    def make_output(self, out_dir, out_prefix):
        """
        Makes full output to the given directory
        """
        out_links = os.path.join(out_dir, out_prefix + "_scaffolds.links")
        out_agp = os.path.join(out_dir, out_prefix + "_scaffolds.agp")
        out_chr = os.path.join(out_dir, out_prefix + "_scaffolds.fasta")
        out_unplaced = os.path.join(out_dir, out_prefix + "_unplaced.fasta")

        self._make_unplaced_fasta()
        write_fasta_dict(self.unplaced_fasta, out_unplaced)

        self._fix_gaps()
        self._make_scaffolds_fasta()
        write_fasta_dict(self.scaffolds_fasta, out_chr)

        self._output_agp(out_agp, out_prefix)
        self._print_statistics()
        output_links(self.scaffolds, out_links)

    def _fix_gaps(self):
        """
        Handles negative gaps, ensures that gap values are
        within deined range
        """
        def get_seq(contig):
            seq_name, seg_start, seg_end = contig.name_with_coords()
            cont_seq = self.fragments_fasta[seq_name][seg_start:seg_end]
            if contig.sign < 0:
                cont_seq = reverse_complement(cont_seq)
            return cont_seq

        def count_ns(cnt_1, cnt_2):
            seq_1, seq_2 = get_seq(cnt_1), get_seq(cnt_2)

            left_ns, right_ns = 0, 0
            for i in range(len(seq_1) - 1, 0, -1):
                if seq_1[i].upper() != "N":
                    break
                left_ns += 1
            for i in range(len(seq_2) - 1):
                if seq_2[i].upper() != "N":
                    break
                right_ns += 1
            return left_ns, right_ns

        for scf in self.scaffolds:
            for cnt_1, cnt_2 in zip(scf.contigs[:-1], scf.contigs[1:]):
                if cnt_1.link.supporting_assembly:
                    cnt_1.trim_right(max(0, -cnt_1.link.gap))
                    cnt_1.link.gap = max(0, cnt_1.link.gap)
                    continue

                left_ns, right_ns = count_ns(cnt_1, cnt_2)

                cnt_1.trim_right(left_ns)
                cnt_2.trim_left(right_ns)
                cnt_1.link.gap += left_ns + right_ns
                cnt_1.link.gap = max(cnt_1.link.gap,
                                     config.vals["min_scaffold_gap"])
                #cnt_1.link.gap = min(cnt_1.link.gap,
                #                     config.vals["max_scaffold_gap"])

    def _output_agp(self, out_agp, assembly_name):
        """
        Output file in NCBI AGP format
        """
        SHIFT = 1
        with open(out_agp, "w") as f:
            f.write("##agp-version  2.0\n")
            f.write("#ASSEMBLY NAME: {0}\n".format(assembly_name))
            f.write("#DESCRIPTION: Pseudochromosome assembly\n")
            f.write("#PROGRAM: Ragout v{0}\n".format(__version__))
            for scf in sorted(self.scaffolds, key=lambda s: s.name):
                chr_pos = 0
                for contig_id, contig in enumerate(scf.contigs):
                    chr_start = chr_pos
                    chr_end = chr_pos + contig.length()
                    chr_pos = chr_end + contig.link.gap
                    cont_name, cont_start, cont_end = contig.name_with_coords()
                    strand = "+" if contig.sign > 0 else "-"
                    support = _support_to_string(contig.link)

                    contig_num = 2 * contig_id + 1
                    gap_num = 2 * contig_id + 2
                    cont_fields = [scf.name, chr_start + SHIFT, chr_end,
                                   contig_num, "W", cont_name, cont_start + SHIFT,
                                   cont_end, strand]
                    f.write("\t".join(map(str, cont_fields)) + "\n")
                    if contig.link.gap > 0:
                        gap_fields = [scf.name, chr_end + SHIFT, chr_pos, gap_num,
                                      "N", contig.link.gap,
                                      "scaffold", "yes", support]
                        f.write("\t".join(map(str, gap_fields)) + "\n")

    def _make_unplaced_fasta(self):
        """
        Creates unplaced (not used in scaffolds) sequences in FASTA format
        """
        used_ranges_by_seq = defaultdict(list)
        for scf in self.scaffolds:
            for ctg in scf.contigs:
                seq_name, seq_start, seq_end = ctg.name_with_coords()
                used_ranges_by_seq[seq_name].append((seq_start, seq_end))
        for seq_name in self.fragments_fasta:
            seq_len = len(self.fragments_fasta[seq_name])
            used_ranges_by_seq[seq_name].append((0, 0))
            used_ranges_by_seq[seq_name].append((seq_len, seq_len))
            used_ranges_by_seq[seq_name].sort()

        unused_ranges_by_seq = defaultdict(list)
        for seq_name in self.fragments_fasta:
            for range_1, range_2 in zip(used_ranges_by_seq[seq_name][:-1],
                                        used_ranges_by_seq[seq_name][1:]):
                if range_1[1] < range_2[0]:
                    unused_ranges_by_seq[seq_name].append((range_1[1],
                                                           range_2[0]))

        unplaced_fasta = {}
        for seq_name, unused_ranges in unused_ranges_by_seq.items():
            for ur in unused_ranges:
                if ur[0] == 0 and ur[1] == len(self.fragments_fasta[seq_name]):
                    fragment_name = seq_name
                else:
                    fragment_name = seq_name + "[{0}:{1}]".format(ur[0], ur[1])
                unplaced_fasta[fragment_name] = \
                                    self.fragments_fasta[seq_name][ur[0]:ur[1]]

        self.unplaced_fasta = unplaced_fasta

    def _make_scaffolds_fasta(self):
        """
        Creates FASTA sequences of scaffolds
        """
        logger.info("Generating FASTA output")
        scaffolds_fasta = {}

        self.used_fragments_len = 0
        self.introduced_gap_len = 0
        for scf in self.scaffolds:
            scf_seqs = []
            for contig in scf.contigs:
                seq_name, seq_start, seq_end = contig.name_with_coords()
                cont_seq = self.fragments_fasta[seq_name][seq_start:seq_end]
                if contig.sign < 0:
                    cont_seq = reverse_complement(cont_seq)

                scf_seqs.append(cont_seq)
                if contig.link.gap > 0:
                    scf_seqs.append("N" * contig.link.gap)
                    self.introduced_gap_len += contig.link.gap

                self.used_fragments_len += len(cont_seq)

            scf_seq = "".join(scf_seqs)
            scaffolds_fasta[scf.name] = scf_seq

        self.scaffolds_fasta = scaffolds_fasta

    def _print_statistics(self):
        """
        Computes and prints some useful statistics
        """
        unplaced_len = sum(map(len, list(self.unplaced_fasta.values())))
        fragments_len = sum(map(len, list(self.fragments_fasta.values())))
        output_len = self.used_fragments_len + self.introduced_gap_len

        #used_perc = 100 * float(self.used_fragments_len) / fragments_len
        unplaced_perc = 100 * float(unplaced_len) / fragments_len
        gap_perc = 100 * float(self.introduced_gap_len) / output_len

        unplaced_count = len(self.unplaced_fasta)
        used_fragments_num = 0
        for scf in self.scaffolds:
            used_fragments_num += len(scf.contigs)

        contigs_len = [len(c) for c in self.fragments_fasta.values()]
        scaffolds_len = [len(c) for c in self.scaffolds_fasta.values()]
        contigs_n50 = _calc_n50(contigs_len, fragments_len)
        scaffolds_n50 = _calc_n50(scaffolds_len, output_len)

        logger.info("Assembly statistics:\n\n"
                    "\tScaffolds:\t\t{0}\n"
                    "\tUsed fragments:\t\t{1}\n"
                    "\tScaffolds length:\t{2}\n\n"
                    "\tUnplaced fragments:\t{3}\n"
                    "\tUnplaced length:\t{4} ({5:2.2f}%)\n"
                    "\tIntroduced Ns length:\t{6} ({7:2.2f}%)\n\n"
                    "\tFragments N50:\t\t{8}\n"
                    "\tAssembly N50:\t\t{9}\n"
                    .format(len(self.scaffolds), used_fragments_num,
                            output_len, unplaced_count, unplaced_len,
                            unplaced_perc, self.introduced_gap_len, gap_perc,
                            contigs_n50, scaffolds_n50))


def output_links(scaffolds, out_links):
    """
    Outputs pretty table with information about adjacencies
    """
    HEADER = ["sequence", "start", "length", "gap", "support"]
    COL_GAP = 4

    with open(out_links, "w") as f:
        for scf in sorted(scaffolds, key=lambda s: s.name):
            rows = []
            cur_pos = 0

            for contig in scf.contigs:
                start = cur_pos
                cur_pos = start + contig.length() + contig.link.gap
                support = _support_to_string(contig.link)

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


def _support_to_string(link):
    """
    Converts information about supporting adjacencies to string.
    Could be used separately form OutputGenerator for debugging purposes
    """
    supp_genomes = sorted(link.supporting_genomes)
    support_to_str = lambda gc: "{0}:{1}".format(gc.genome, gc.chr)
    support = ",".join(map(support_to_str, supp_genomes))
    if link.supporting_assembly:
        support += ",~>"
    return support


def _calc_n50(scaffolds_lengths, assembly_len):
    n50 = 0
    sum_len = 0
    for l in sorted(scaffolds_lengths, reverse=True):
        sum_len += l
        if sum_len > assembly_len // 2:
            n50 = l
            break
    return n50
