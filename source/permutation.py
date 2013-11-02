from collections import namedtuple, defaultdict
import os


class Permutation:
    def __init__(self, ref_id, chr_id, chr_num, blocks):
        self.ref_id = ref_id
        self.chr_id = chr_id
        self.chr_num = chr_num
        self.blocks = blocks

    def iter_blocks(self, circular=False):
        assert len(self.blocks) > 0
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


#TODO: place somewhere else
def parse_config(filename):
    prefix = os.path.dirname(filename)
    references = {}
    target = {}
    tree_str = None
    block_size = None

    for line in open(filename, "r").read().splitlines():
        if line.startswith("#"):
            continue

        if line.startswith("REF"):
            ref_id, ref_file = line[4:].split("=")
            references[ref_id] = os.path.join(prefix, ref_file)

        if line.startswith("TARGET"):
            ref_id, ref_file = line[7:].split("=")
            target[ref_id] = os.path.join(prefix, ref_file)

        if line.startswith("TREE"):
            tree_str = line.split("=")[1]

        if line.startswith("BLOCK"):
            sizes = line.split("=")[1].split(",")
            block_size = map(int, sizes)

    return references, target, tree_str, block_size


class PermutationContainer:
    def __init__(self, config_file):
        self.ref_perms = []
        self.target_perms = []

        ref_files, target_files, _tree, _blocks = parse_config(config_file)
        for ref_id, ref_file in ref_files.iteritems():
            self.ref_perms.extend(parse_blocks_file(ref_id, ref_file))

        for t_id, t_file in target_files.iteritems():
            self.target_perms.extend(parse_blocks_file(t_id, t_file))

        self.duplications = find_duplications(self.ref_perms, self.target_perms)

        self.target_blocks = set()
        for perm in self.target_perms:
            self.target_blocks |= set(map(abs, perm.blocks))

        to_hold = self.target_blocks - self.duplications
        self.ref_perms_filtered = [filter_perm(p, to_hold) for p in self.ref_perms]
        self.target_perms_filtered = [filter_perm(p, to_hold) for p in self.target_perms]
