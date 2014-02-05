from collections import defaultdict
from permutations import Block

def maf_to_permutations(maf_file, min_block):
    filename = maf_file
    block_id = 0
    permutations = defaultdict(list)
    new_lcb = False

    for line in open(filename, "r").read().splitlines():
        if line.startswith("a"):
            new_lcb = True
        elif line.startswith("s"):
            vals = line.split("\t")
            genome = vals[1]
            start = int(vals[2])
            length = int(vals[3])
            strand = vals[4]
            src_len = int(vals[5])
            absolute_start = start if strand == "+" else src_len - (start + length)

            if length > min_block:
                if new_lcb:
                    block_id += 1
                    new_lcb = False
                numeric_id = int("{0}{1}".format(strand, block_id))
                permutations[genome].append(Block(numeric_id, absolute_start, length))

    for lst in permutations.values():
        lst.sort(key=lambda b: b.start)
    return permutations
