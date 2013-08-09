import os
from collections import namedtuple, defaultdict

SyntenyBlock = namedtuple("SyntenyBlock", ["seq", "chr_id", "strand", "id", "start", "end", "chr_num"])
SeqInfo = namedtuple("SeqInfo", ["id", "length"])
Permutation = namedtuple("Permutation", ["chr_id", "chr_num", "blocks"])

class Contig:
    def __init__(self, name):
        self.name = name
        self.sign = 1
        self.blocks = []


class SibeliaOutput:
    def __init__(self, sibelia_dir, contig_names):
        self.parse_sibelia_output(sibelia_dir, contig_names)

    def parse_sibelia_output(self, sibelia_dir, contig_names):
        coords_file = os.path.join(sibelia_dir, "blocks_coords.txt")
        permutations_file = os.path.join(sibelia_dir, "genomes_permutations.txt")

        blocks_info, seq_info = self.parse_coords_file(coords_file)
        self.blocks_info = blocks_info
        self.seq_info = seq_info

        permutations, contigs = self.parse_permutations_file(permutations_file, contig_names)
        self.permutations = permutations
        self.contigs = contigs
        #return permutations, contigs, blocks_coords


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


    def get_reference_count(self):
        return len(self.permutations)


    def build_contig_index(self):
        index = defaultdict(list)
        for i, c in enumerate(self.contigs):
            for block in c.blocks:
                index[abs(block)].append(i)
        return index


    def get_blocks_distance(self, left_block, right_block, ref_num):
        """
        Only non-duplicated blocks
        """
        blocks_coords = self.blocks_info
        left_instances = filter(lambda b: b.chr_num == ref_num, blocks_coords[left_block])
        right_instances = filter(lambda b: b.chr_num == ref_num, blocks_coords[right_block])

        assert len(left_instances) == len(right_instances) == 1
        if left_instances[0].strand == "+":
            left = left_instances[0].end
        else:
            left = left_instances[0].start

        if right_instances[0].strand == "+":
            right = right_instances[0].start
        else:
            right = right_instances[0].end

        assert right >= left
        return right - left - 1


    def block_offset(self, block, contig_name, reverse=False):
        """
        distance from beginning of contig to block, or from end of contig to block,
        if reverse is specified.
        non-duplicated blocks
        """
        blocks_coords = self.blocks_info
        instances = filter(lambda b: b.chr_id == contig_name, blocks_coords[block])
        #print blocks_coords[block]
        #print block, contig_name
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
        #print "-", reverse, block, offset
        return offset

