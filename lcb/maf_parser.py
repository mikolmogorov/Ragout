from collections import defaultdict
from permutations import Block

def maf_to_permutations(maf_file, min_block):
    MAX_GAP_RATE = 0.3

    block_id = 0
    permutations = defaultdict(list)
    new_blocks = defaultdict(list)
    new_lcb = False

    def update_perms():
        if len(new_blocks) > 1:
            for seq_id in new_blocks:
                permutations[seq_id].extend(new_blocks[seq_id])
        new_blocks.clear()

    for line in open(maf_file, "r").read().splitlines():
        if line.startswith("a"):
            new_lcb = True
            update_perms()

        elif line.startswith("s"):
            seq_id, start, ungapped_len, strand, src_len, seq = line.split("\t")[1:]
            start = int(start)
            src_len = int(src_len)
            ungapped_len = int(ungapped_len)
            gapped_len = len(seq)
            absolute_start = start if strand == "+" else src_len - (start + ungapped_len)

            gap_rate = float(gapped_len - ungapped_len) / float(gapped_len)

            if ungapped_len > min_block and gap_rate < MAX_GAP_RATE:
                if new_lcb:
                    block_id += 1
                    new_lcb = False
                numeric_id = int("{0}{1}".format(strand, block_id))
                new_blocks[seq_id].append(Block(numeric_id, absolute_start, ungapped_len))

    update_perms()

    for lst in permutations.values():
        lst.sort(key=lambda b: b.start)
    return permutations
