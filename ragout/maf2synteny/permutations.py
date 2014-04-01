from __future__ import print_function
from collections import namedtuple, defaultdict
from bisect import bisect_right
from copy import copy

Block = namedtuple("Block", ["id", "start", "length"])

def output_permutations(permutations, stream):
    for seq_id, blocks in permutations.iteritems():
        stream.write(">{0}\n".format(seq_id))
        for block in blocks:
            sign = "+" if block.id > 0 else "-"
            stream.write("{0}{1} ".format(sign, abs(block.id)))
        stream.write("$\n")


def output_blocks_coords(permutations, seq_length, stream):
    num_ids = dict(map(reversed, enumerate(permutations.keys())))
    by_block = defaultdict(list)

    for seq_id, blocks in permutations.iteritems():
        for block in blocks:
            by_block[abs(block.id)].append((block, num_ids[seq_id]))

    stream.write("Seq_id\tSize\tDescription\n")
    for seq, id in num_ids.iteritems():
        stream.write("{0}\t{1}\t{2}\n".format(id, seq_length[seq], seq))
    stream.write("-" * 80 + "\n")

    for block_id, blocklist in by_block.iteritems():
        blocklist.sort(key=lambda b: b[1])
        stream.write("Block #{0}\nSeq_id\tStrand\tStart\tEnd\tLength\n"
                     .format(block_id))

        for block, seq_id in blocklist:
            strand = "+" if block.id > 0 else "-"
            start = block.start
            length = block.length
            end = start + length
            stream.write("{0}\t{1}\t{2}\t{3}\t{4}\n"
                         .format(seq_id, strand, start, end, length))
        stream.write("-" * 80 + "\n")


def load_permutations(file):
    permutations = {}
    for line in open(file, "r").read().splitlines():
        if line.startswith(">"):
            seq_id = line[1:]
        else:
            blocks = map(lambda i: Block(int(i), 0, 0), line.split(" ")[:-1])
            permutations[seq_id] = blocks

    return permutations


def renumerate(permutations):
    next_id = [1]
    new_ids = {}
    def new_id(old_id):
        sign = 1 if old_id > 0 else -1
        if abs(old_id) not in new_ids:
            new_ids[abs(old_id)] = next_id[0]
            next_id[0] += 1
        return sign * new_ids[abs(old_id)]

    permutations = copy(permutations)
    new_block = lambda b: Block(new_id(b.id), b.start, b.length)
    for seq_id in permutations:
        permutations[seq_id] = map(new_block, permutations[seq_id])
    return permutations


def filter_by_size(permutations, min_size, min_flank=0, block_groups={}):
    permutations = copy(permutations)
    block_to_group = defaultdict(lambda: None)
    for group, blocks in block_groups.iteritems():
        for block in blocks:
            block_to_group[block] = group
    group_len = defaultdict(int)

    for seq_id, blocks in permutations.iteritems():
        for block in blocks:
            group_id = block_to_group[abs(block.id)]
            if not group_id:
                continue
            group_len[(group_id, seq_id)] += block.length

    should_output = set()
    for seq_id, blocks in permutations.iteritems():
        for block in blocks:
            group_id = block_to_group[abs(block.id)]
            if block.length >= min_size:
                should_output.add(abs(block.id))
            else:
                group_id = block_to_group[abs(block.id)]
                if (group_len[(group_id, seq_id)] >= min_size and
                        block.length > min_flank):
                    should_output.add(abs(block.id))

    for gen in permutations.keys():
        new_blocks = [b for b in permutations[gen] if abs(b.id) in should_output]
        if len(new_blocks):
            permutations[gen] = new_blocks
        else:
            del permutations[gen]

    return permutations


def output_statistics(permutations, seq_length, stream):
    multiplicity = defaultdict(int)
    num_ids = dict(map(reversed, enumerate(permutations.keys())))
    by_block = defaultdict(list)
    covered = defaultdict(int)

    for seq_id, blocks in permutations.iteritems():
        for block in blocks:
            covered[seq_id] += block.length
            by_block[abs(block.id)].append((block, num_ids[seq_id]))

    stream.write("Seq_id\tSize\tDescription\n")
    for seq, id in num_ids.iteritems():
        stream.write("{0}\t{1}\t{2}\n".format(id, seq_length[seq], seq))
    stream.write("-" * 80 + "\n")

    for block_id, blocklist in by_block.iteritems():
        multiplicity[len(blocklist)] += 1

    for mul, num in sorted(multiplicity.iteritems()):
        stream.write("{0}\t{1}\n".format(mul, num))
    stream.write("-" * 80 + "\n")

    for seq_id, cov in covered.iteritems():
        stream.write("{0}\t{1:4.2f}%\n"
                      .format(seq_id, 100 * float(cov) / seq_length[seq_id]))


def merge_permutations(simplified_perms, initial_perms):
    permutations = defaultdict(list, simplified_perms)
    by_block = defaultdict(list)
    starts, ends = defaultdict(list), defaultdict(list)
    next_id = 0

    for seq_id, blocks in simplified_perms.iteritems():
        starts[seq_id] = map(lambda b: b.start, blocks)
        ends[seq_id] = map(lambda b: b.start + b.length, blocks)

    for seq_id, blocks in initial_perms.iteritems():
        for block in blocks:
            by_block[abs(block.id)].append((block, seq_id))
            next_id = max(next_id, abs(block.id))
    next_id += 1

    for block_id, blocklist in by_block.iteritems():
        inserted_blocks = []
        for block, seq_id in blocklist:
            left_i = bisect_right(ends[seq_id], block.start)
            right_i = bisect_right(starts[seq_id], block.start + block.length)
            if left_i == right_i:
                inserted_blocks.append((block, seq_id))
                #if 0 < left_i < len(starts[seq_id]):
                    #x1 = (simplified_perms[seq_id][left_i - 1].start +
                    #      simplified_perms[seq_id][left_i - 1].length)
                    #x2, x3 = block.start, block.start + block.length
                    #x4 = simplified_perms[seq_id][left_i].start
                    #assert x1 <= x2 and x3 <= x4

        if len(inserted_blocks) > 1:
            for block, seq_id in inserted_blocks:
                sign = 1 if block.id > 0 else -1
                permutations[seq_id].append(Block(next_id * sign, block.start,
                                                  block.length))
            next_id += 1

    for seq_id in permutations:
        permutations[seq_id].sort(key=lambda b: b.start)

    return permutations
