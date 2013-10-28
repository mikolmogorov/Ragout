import os
import shutil
import subprocess
import copy
from collections import namedtuple, defaultdict
from datatypes import *
from Bio import SeqIO
from Bio.Seq import Seq


SyntenyBlock = namedtuple("SyntenyBlock", ["seq", "chr_id", "strand", "id", "start", "end", "chr_num"])
SeqInfo = namedtuple("SeqInfo", ["id", "length"])
Permutation = namedtuple("Permutation", ["chr_id", "chr_num", "blocks"])


SIBELIA_BIN = "/home/volrath/Bioinf/Sibelia/distr/bin/Sibelia"


def run_sibelia(references, contigs, block_size, out_dir):
    print "Running Sibelia..."
    cmdline = [SIBELIA_BIN, "-s", "loose", "-m", str(block_size), "-o", out_dir]
    cmdline.extend(references)
    cmdline.append(contigs)
    process = subprocess.Popen(cmdline)
    process.wait()
    os.remove(os.path.join(out_dir, "coverage_report.txt"))
    os.remove(os.path.join(out_dir, "d3_blocks_diagram.html"))
    shutil.rmtree(os.path.join(out_dir, "circos"))


def build_contig_index(contigs):
    index = defaultdict(list)
    for i, c in enumerate(contigs):
        for block in c.blocks:
            index[abs(block)].append(i)
    return index


class SibeliaRunner:
    def __init__(self, references, target_file, block_size, output_dir, contig_names, skip_run):
        if not skip_run:
            run_sibelia(references.values(), target_file, block_size, output_dir)
        self.parse_references(references)
        self.parse_sibelia_output(output_dir, contig_names, references)


    def parse_sibelia_output(self, sibelia_dir, contig_names, references):
        coords_file = os.path.join(sibelia_dir, "blocks_coords.txt")
        permutations_file = os.path.join(sibelia_dir, "genomes_permutations.txt")

        blocks_info, seq_info = self.parse_coords_file(coords_file)
        self.blocks_info = blocks_info
        self.seq_info = seq_info

        permutations, contigs = self.parse_permutations_file(permutations_file, contig_names)
        self.permutations = permutations
        self.contigs = contigs
        self.references = references


    def parse_references(self, references):
        chr_id_to_ref_id = {}
        for ref_id, filename in references.iteritems():
            seq = SeqIO.parse(filename, "fasta")
            chr_id_to_ref_id[seq.next().id] = ref_id
        self.chr_id_to_ref_id = chr_id_to_ref_id


    def get_filtered_contigs(self):
        dups = self.get_duplications()

        new_contigs = []
        for contig in self.contigs:
            new_blocks = [b for b in contig.blocks if abs(b) not in dups]
            if new_blocks:
                new_contigs.append(copy.copy(contig))
                new_contigs[-1].blocks = new_blocks

        return new_contigs


    def get_duplications(self):
        duplications = set()
        for perm in self.permutations:
            current = set()
            for block in perm.blocks:
                if abs(block) in current:
                    duplications.add(abs(block))
                current.add(abs(block))

        current = set()
        for contig in self.contigs:
            for block in contig.blocks:
                if abs(block) in current:
                    duplications.add(abs(block))
                current.add(abs(block))

        return duplications


    def parse_permutations_file(self, filename, contig_names):
        fin = open(filename, "r")
        contigs = []
        permutations = []
        contig_name = None
        ref_name = None
        num_by_id = {seq.id : seq_num for seq_num, seq in self.seq_info.iteritems()}

        for line in fin:
            if line.startswith(">"):
                name = line.strip()[1:]
                if name in contig_names:
                    contig_name = name
                else:
                    ref_name = name
                continue

            blocks = line.strip().split(" ")[0:-1]

            #contig
            if contig_name:
                contig = Contig(contig_name)
                contig.blocks = map(int, blocks)
                contigs.append(contig)
            #reference
            else:
                ref_num = num_by_id[ref_name]
                permutations.append(Permutation(chr_id=ref_name, chr_num=ref_num,
                                                blocks=map(int, blocks)))
        return (permutations, contigs)


    def parse_coords_file(self, blocks_file):
        group = [[]]
        seq_info = {}
        blocks_info = {}
        line = [l.strip() for l in open(blocks_file) if l.strip()]
        for l in line:
            if l.startswith("-"):
                group.append([])
            else:
                group[-1].append(l)
        for l in group[0][1:]:
            seq_num, seq_len, seq_id = l.split()
            seq_num = int(seq_num)
            seq_info[seq_num] = SeqInfo(seq_id, int(seq_len))
        for g in [g for g in group[1:] if g]:
            block_id = int(g[0].split()[1][1:])
            blocks_info[block_id] = []
            for l in g[2:]:
                chr_num, bl_strand, bl_start, bl_end, bl_length = l.split()
                chr_num = int(chr_num)
                chr_id = seq_info[chr_num].id
                blocks_info[block_id].append(SyntenyBlock(seq="", chr_id=chr_id,
                                            strand=bl_strand, id=block_id,
                                            start=int(bl_start), end=int(bl_end),
                                            chr_num=int(chr_num)))
        return (blocks_info, seq_info)


    #def get_references_count(self):
    #    return len(self.permutations)

"""
    def get_blocks_distance(self, left_block, right_block, ref_num):
        blocks_coords = self.blocks_info
        left_instances = filter(lambda b: b.chr_num == ref_num, blocks_coords[left_block])
        right_instances = filter(lambda b: b.chr_num == ref_num, blocks_coords[right_block])

        #Only non-duplicated blocks
        assert len(left_instances) == len(right_instances) == 1
        if left_instances[0].strand == "+":
            left = left_instances[0].end
        else:
            left = left_instances[0].start

        if right_instances[0].strand == "+":
            right = right_instances[0].start
        else:
            right = right_instances[0].end

        #assert right >= left
        #TODO: correct distance here
        return abs(right - left) - 1


    def block_offset(self, block, contig_name, reverse=False):
        #distance from beginning of contig to block, or from end of contig to block,
        #if reverse is specified.
        #non-duplicated blocks
        blocks_coords = self.blocks_info
        instances = filter(lambda b: b.chr_id == contig_name, blocks_coords[block])
        assert len(instances) == 1

        contig_len = self.seq_info[instances[0].chr_num].length
        if not reverse:
            if instances[0].strand == "+":
                offset = instances[0].start
            else:
                offset = instances[0].end
        else:
            if instances[0].strand == "+":
                offset = contig_len - instances[0].end
            else:
                offset = contig_len - instances[0].start
        return offset
"""
