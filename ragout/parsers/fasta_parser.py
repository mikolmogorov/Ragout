from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def read_fasta_dict(filename):
    fasta_dict = {}
    for seq in SeqIO.parse(filename, "fasta"):
        fasta_dict[seq.id] = str(seq.seq)
    return fasta_dict


def write_fasta_dict(fasta_dict, filename):
    out = open(filename, "w")
    for name, seq in fasta_dict.items():
        SeqIO.write(SeqRecord(Seq(seq), id=name, description=""),
                    out, "fasta")

def reverse_complement(string):
    return str(Seq(string).reverse_complement())
