from collections import namedtuple, defaultdict

AlignmentInfo = namedtuple("AlignmentInfo", ["s_ref", "e_ref", "s_qry", "e_qry",
                            "len_ref", "len_qry", "ref_id", "contig_id"])

class Hit:
    def __init__(self, index, chr, coord, sign):
        self.index = index
        self.chr = chr
        self.coord = coord
        self.sign = sign

    def __str__(self):
        return str(self.index) + " : " + str(self.chr)


def group_by_chr(alignment):
    by_chr = defaultdict(list)
    for entry in alignment:
        by_chr[entry.ref_id].append(entry)
    for chr_id in by_chr:
        by_chr[chr_id].sort(key=lambda e: e.s_ref)
    return by_chr


def get_order(alignment):
    chr_len = defaultdict(int)
    contig_len = defaultdict(int)

    for entry in alignment:
        contig_len[entry.contig_id] = max(entry.len_qry, contig_len[entry.contig_id])
        chr_len[entry.ref_id] = max(entry.s_ref, chr_len[entry.ref_id])

    by_chr = group_by_chr(alignment)
    entry_ord = defaultdict(list)

    contig_pos = 1
    prev_start = None
    for chr_id, alignment in by_chr.items():
        for e in alignment:
            if prev_start is not None and e.s_ref > prev_start:
                    contig_pos += 1
            prev_start = e.s_ref
            sign = 1 if e.e_qry > e.s_qry else -1
            entry_ord[e.contig_id].append(Hit(contig_pos, chr_id, e.s_ref, sign))

    return entry_ord, chr_len, contig_len


def join_collinear(alignment):
    new_entries = []
    by_chr = group_by_chr(alignment)
    for chr_id in by_chr:
        by_chr[chr_id].sort(key=lambda e: e.s_ref)
        #prev_contig = None
        start_entry = None
        last_entry = None
        for entry in by_chr[chr_id]:
            if not start_entry:
                start_entry = entry
            elif start_entry.contig_id != entry.contig_id:
                new_entries.append(AlignmentInfo(start_entry.s_ref, last_entry.e_ref,
                                    start_entry.s_qry, last_entry.e_qry,
                                    abs(last_entry.e_ref - start_entry.s_ref),
                                    abs(last_entry.e_qry - start_entry.s_qry),
                                    last_entry.ref_id, last_entry.contig_id))
                start_entry = entry
            last_entry = entry

    return new_entries


def parse_nucmer_coords(filename):
    chr_alias = {}
    chr_num = 1

    alignment = []
    for line in open(filename, "r"):
        line = line.strip()
        if not len(line) or not line[0].isdigit():
            continue

        vals = line.split(" | ")
        s_ref, e_ref = list(map(int, vals[0].split()))
        s_qry, e_qry = list(map(int, vals[1].split()))
        len_ref, len_qry = list(map(int, vals[2].split()))
        ref_id, contig_id = vals[4].split("\t")

        if ref_id not in chr_alias:
            chr_alias[ref_id] = "chr{0}".format(chr_num)
            chr_num += 1
        alignment.append(AlignmentInfo(s_ref, e_ref, s_qry, e_qry,
                            len_ref, len_qry, chr_alias[ref_id], contig_id))

    return alignment


def filter_by_coverage(alignment):
    MIN_HIT = 0.45
    by_name = defaultdict(list)
    for entry in alignment:
        by_name[entry.contig_id].append(entry)

    for name in by_name:
        by_name[name].sort(key=lambda e: e.len_qry, reverse=True)
        len_filter = lambda e: e.len_qry > MIN_HIT * by_name[name][0].len_qry
        by_name[name] = list(filter(len_filter, by_name[name]))
        #print(by_name[name])

    filtered_alignment = []
    for ent_lst in by_name.values():
        filtered_alignment.extend(ent_lst)

    return filtered_alignment
