from collections import namedtuple, defaultdict

SyntenyBlock = namedtuple("SyntenyBlock", ["seq", "chr_id", "strand", "id", "start", "end", "chr_num"])
Permutation = namedtuple("Permutation", ["chr_id", "chr_num", "blocks"])

class Contig:
    def __init__(self, name):
        self.name = name
        self.sign = 1
        self.blocks = []


def parse_permutations_file(filename):
    fin = open(filename, "r")
    contigs = []
    permutations = []
    contig_name = None
    ref_name = None
    ref_num = 0

    for line in fin:
        if line.startswith(">"):
            if line.startswith(">contig") or line.startswith(">scaffold"):
                contig_name = line.strip()[1:]
            else:
                ref_name = line.strip()[1:]
            continue

        blocks = line.strip().split(" ")[0:-1]

        #contig
        if contig_name:
            contig = Contig(contig_name)
            contig.blocks = map(int, blocks)
            contigs.append(contig)
        #reference
        else:
            permutations.append(Permutation(chr_id=ref_name, chr_num=ref_num,
                                            blocks=map(int, blocks)))
    return (permutations, contigs)


def parse_coords_file(blocks_file):
    group = [[]]
    num_seq_id = dict()
    seq_id_num = dict()
    line = [l.strip() for l in open(blocks_file) if l.strip()]
    for l in line:
        if l[0] == "-":
            group.append([])
        else:
            group[-1].append(l)
    for l in group[0][1:]:
        l = l.split()
        num_seq_id[l[0]] = l[2]
        seq_id_num[l[2]] = int(l[0])
    ret = dict()
    for g in [g for g in group[1:] if g]:
        block_id = int(g[0].split()[1][1:])
        ret[block_id] = []
        for l in g[2:]:
            l = l.split()
            chr_id = num_seq_id[l[0]]
            start = int(l[2])
            end = int(l[3])
            chr_num = int(l[0])
            ret[block_id].append(SyntenyBlock(seq="", chr_id=chr_id, strand=l[1],
                                id=block_id, start=start, end=end, chr_num=(chr_num - 1)))
    return (ret, seq_id_num)


def build_contig_index(contigs):
    index = defaultdict(list)
    for i, c in enumerate(contigs):
        for block in c.blocks:
            index[abs(block)].append(i)
    return index


def get_blocks_distance(left_block, right_block, ref_num, blocks_coords):
    """
    Only non-duplicated blocks
    """
    left_instances = filter(lambda b: b.chr_num == ref_num, blocks_coords[left_block])
    right_instances = filter(lambda b: b.chr_num == ref_num, blocks_coords[right_block])

    #print len(left_instances), len(right_instances), right_instances
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

