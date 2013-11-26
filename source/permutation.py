from collections import defaultdict
import config_parser as parser


class Permutation:
    def __init__(self, ref_id, chr_id, chr_num, blocks):
        self.ref_id = ref_id
        self.chr_id = chr_id
        self.chr_num = chr_num
        self.blocks = blocks
        self.target_perms = []
        self.ref_perms = []
        self.ref_perms_filtered = []
        self.target_perms_filtered = []

    def iter_blocks(self, circular=False):
        if not len(self.blocks):
            return

        for block in self.blocks:
            yield block

        if circular:
            yield self.blocks[0]


def find_duplications(ref_perms, target_perms):
    index = defaultdict(set)
    duplications = set()
    for perm in ref_perms + target_perms:
        for block in map(abs, perm.blocks):
            if perm.ref_id in index[block]:
                duplications.add(block)
            else:
                index[block].add(perm.ref_id)

    return duplications


def filter_perm(perm, to_hold):
    new_perm = Permutation(perm.ref_id, perm.chr_id, perm.chr_num, [])
    for block in perm.blocks:
        if abs(block) in to_hold:
            new_perm.blocks.append(block)
    return new_perm


def parse_blocks_file(ref_id, filename):
    name = ""
    permutations = []
    chr_count = 0
    for line in open(filename, "r").read().splitlines():
        if line.startswith(">"):
            name = line[1:]
        else:
            blocks = line.split(" ")[:-1]
            permutations.append(Permutation(ref_id, name, chr_count, map(int, blocks)))
            chr_count += 1
    return permutations


class PermutationContainer:
    def __init__(self, config_file):
        self.ref_perms = []
        self.target_perms = []

        config = parser.parse_ragout_config(config_file)
        for ref_id, ref_file in config.references.iteritems():
            self.ref_perms.extend(parse_blocks_file(ref_id, ref_file))

        for t_id, t_file in config.targets.iteritems():
            self.target_perms.extend(parse_blocks_file(t_id, t_file))

        self.duplications = find_duplications(self.ref_perms, self.target_perms)

        self.target_blocks = set()
        for perm in self.target_perms:
            self.target_blocks |= set(map(abs, perm.blocks))

        to_hold = self.target_blocks - self.duplications
        self.ref_perms_filtered = [filter_perm(p, to_hold) for p in self.ref_perms]
        self.target_perms_filtered = [filter_perm(p, to_hold) for p in self.target_perms]
