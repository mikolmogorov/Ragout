from collections import defaultdict, namedtuple
from permutations import Block
from itertools import izip
import sys

Hit = namedtuple("Column", ["seq_id", "start", "seq_len", "len", "strand", "seq"])


def extend_column(prev_col, next_col, max_ref_gap):
    if len(prev_col) != len(next_col):
        return None

    new_column = []
    for num, (prev_hit, next_hit) in enumerate(izip(prev_col, next_col)):
        if prev_hit.seq_id != next_hit.seq_id:
            return None
        if prev_hit.strand != next_hit.strand:
            return None

        diff = next_hit.start - (prev_hit.start + prev_hit.len)
        #num == 0 -> reference
        if ((diff != 0 and num != 0) or
            (num == 0 and not (0 <= diff <= max_ref_gap))):
            return None

        new_len = next_hit.len + prev_hit.len + diff
        new_seq = prev_hit.seq + "-" * diff + next_hit.seq
        new_column.append(Hit(prev_hit.seq_id, prev_hit.start, prev_hit.seq_len,
                              new_len, prev_hit.strand, new_seq))

    return new_column


def output_column(column, stream):
    assert column
    stream.write("a\n")
    for hit in column:
        stream.write("s\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".
            format(hit.seq_id, hit.start, hit.len,
                   hit.strand, hit.seq_len, hit.seq))

    stream.write("\n")


def condense_maf(input_maf, output_maf, max_ref_gap):
    prev_column = []
    current_column = []
    output_maf = open(output_maf, "w")

    for line in open(input_maf, "r"):
        line = line.strip()
        if line.startswith("a"):
            if not current_column:
                continue

            if prev_column:
                result = extend_column(prev_column, current_column, max_ref_gap)
            else:
                result = current_column

            if not result:  #merge failed
                output_column(prev_column, output_maf)
                prev_column = current_column
            else:
                prev_column = result
            current_column = []

        elif line.startswith("s"):
            seq_id, start, ungap_len, strand, seq_len, seq = line.split("\t")[1:]
            start = int(start)
            seq_len = int(seq_len)
            ungap_len = int(ungap_len)
            gapped_len = len(seq)
            absolute_start = start if strand == "+" else seq_len - (start + ungap_len)
            current_column.append(Hit(seq_id, start, seq_len,
                                        ungap_len, strand, seq))

        elif line:
            output_maf.write(line + "\n")

    result = extend_column(prev_column, current_column, 5)
    if not result:
        output_column(prev_column, output_maf)
        output_column(current_column, output_maf)
    else:
        output_column(result, output_maf)



def maf_to_permutations(maf_file, min_block):
    MAX_GAP_RATE = 0.3

    block_id = 0
    permutations = defaultdict(list)
    new_blocks = defaultdict(list)
    seq_length = {}
    new_lcb = False

    def update_perms():
        if len(new_blocks) > 1:
            for seq_id in new_blocks:
                permutations[seq_id].extend(new_blocks[seq_id])
        new_blocks.clear()

    for line in open(maf_file, "r"):
        line = line.strip()
        if line.startswith("a"):
            new_lcb = True
            update_perms()

        elif line.startswith("s"):
            seq_id, start, ungapped_len, strand, src_len, seq = line.split("\t")[1:]
            if seq_id.startswith("gi") and seq_id[-1] != "|": #fix for cactus
                seq_id += "|"
            start = int(start)
            src_len = int(src_len)
            ungapped_len = int(ungapped_len)
            gapped_len = len(seq)
            absolute_start = start if strand == "+" else src_len - (start + ungapped_len)

            gap_rate = float(gapped_len - ungapped_len) / float(gapped_len)
            if ungapped_len >= min_block and gap_rate < MAX_GAP_RATE:
                if new_lcb:
                    block_id += 1
                    new_lcb = False
                numeric_id = int("{0}{1}".format(strand, block_id))
                new_blocks[seq_id].append(Block(numeric_id, absolute_start, ungapped_len))

            seq_length[seq_id] = src_len

    update_perms()

    for lst in permutations.values():
        lst.sort(key=lambda b: b.start)
    return permutations, seq_length
