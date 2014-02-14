from collections import namedtuple, defaultdict

Block = namedtuple("Block", ["id", "start", "length"])

def output_permutations(permutations, stream):
    for seq_id, blocks in permutations.iteritems():
        stream.write(">{0}\n".format(seq_id))
        for block in blocks:
            sign = "+" if block.id > 0 else "-"
            stream.write("{0}{1} ".format(sign, abs(block.id)))
        stream.write("$\n")


def output_blocks_coords(permutations, stream):
    #multiplicity = defaultdict(int)

    num_ids = dict(map(reversed, enumerate(permutations.keys())))
    by_block = defaultdict(list)

    stream.write("Seq_id\tSize\tDescription\n")
    for seq, id in num_ids.iteritems():
        stream.write("{0}\t{1}\t{2}\n".format(id, 0, seq))
    stream.write("-" * 80 + "\n")

    for seq_id, blocks in permutations.iteritems():
        for block in blocks:
            by_block[abs(block.id)].append((block, num_ids[seq_id]))

    for block_id, blocklist in by_block.iteritems():
        #multiplicity[len(blocklist)] += 1
        blocklist.sort(key=lambda b: b[1])
        stream.write("Block #{0}\nSeq_id\tStrand\tStart\tEnd\tLength\n".format(block_id))

        for block, seq_id in blocklist:
            strand = "+" if block.id > 0 else "-"
            start = block.start
            length = block.length
            end = start + length
            stream.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(seq_id, strand, start, end, length))
        stream.write("-" * 80 + "\n")

    #for mul, num in sorted(multiplicity.iteritems()):
    #    print mul, num


def load_permutations(file):
    permutations = {}
    for line in open(file, "r").read().splitlines():
        if line.startswith(">"):
            seq_id = line[1:]
        else:
            blocks = map(lambda i: Block(int(i), 0, 0), line.split(" ")[:-1])
            permutations[seq_id] = blocks

    return permutations


def filter_by_size(permutations, min_size):
    should_output = set()
    for blocks in permutations.itervalues():
        for block in blocks:
            if block.length >= min_size:
                should_output.add(abs(block.id))

    for gen in permutations.keys():
        new_blocks = [b for b in permutations[gen] if abs(b.id) in should_output]
        if len(new_blocks):
            permutations[gen] = new_blocks
        else:
            del permutations[gen]
