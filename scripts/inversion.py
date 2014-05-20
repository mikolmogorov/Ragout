import sys
import random
from Bio import SeqIO
from Bio.Seq import Seq

from utils.nucmer_parser import *

def parse_contigs(filename):
    contigs = set()
    for line in open(filename, "r"):
        if not line.startswith(">"):
            name = line.strip("\n").replace("=", "_") #fix for nucmer
            contigs.add(name[1:])
    return contigs

def get_unique_contigs_in_numcher(alignment):
    first_filtration = {}
    for e in alignment:
        if e.contig_id in first_filtration:
            first_filtration[e.contig_id] += 1
        else:
            first_filtration[e.contig_id] = 1

    result = []
    for cont, counts in first_filtration.iteritems():
        if counts == 1:
            result.append(cont)

    return result

def get_contigs_with_length(input, dataset_dict, contigs, treshold):
    result = []
    for name in input:
        if name in contigs:
            start, end  = dataset_dict[name]
            if abs(start - end) >= treshold:
                result.append(name)
    return result

def do_job(nucmer_coords, scaffold, number_of_inv, reference, output_reference, treshold):
    contigs = parse_contigs(scaffold)
    alignment = parse_nucmer_coords(nucmer_coords)
    alignment = join_collinear(alignment)
    alignment = filter_by_coverage(alignment)

    unique_seq = get_unique_contigs_in_numcher(alignment)

    dataset_dict = {}
    for e in alignment:
        dataset_dict[e.contig_id] = (int(e.s_ref), int(e.e_ref))

    unique_seq = get_contigs_with_length(unique_seq, dataset_dict, contigs, treshold)

    for seq in SeqIO.parse(reference, "fasta"):
        refer_name = seq.name
        refer = str(seq.seq)

    for _ in range(number_of_inv):
        contig = random.choice(unique_seq)
        start, end = dataset_dict[contig]
        compl = Seq(refer[(end - 1):(start - 1):-1]).reverse_complement()
        refer = refer[:start:1] + str(compl)[::-1] + refer[end:]

        print(contig + " (" + str(start) + ", " + str(end) +  ")")

        unique_seq.remove(contig)
        if not unique_seq:
            break

    with open(output_reference, "w") as out:
        out.write(">" + refer_name + '\n')
        for i in range(0, len(refer), 50):
            out.write(refer[i:i + 50] + '\n')

def main():
    if len(sys.argv) < 7:
        print("Usage: inversion.py <nucmer_coords> <scaffold> <reference> <output_reference> <number_of_inversion> <treshold>")
        return

    nucmer_coords = sys.argv[1]
    scaffold = sys.argv[2]
    reference = sys.argv[3]
    output_reference = sys.argv[4]
    number_of_inv = int(sys.argv[5])
    treshold = int(sys.argv[6])

    do_job(nucmer_coords, scaffold, number_of_inv, reference, output_reference, treshold)

if __name__ == "__main__":
    main()
