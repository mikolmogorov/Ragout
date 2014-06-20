"""
This module provides some basic FASTA I/O
"""

class FastaError(Exception):
    pass

def read_fasta_dict(filename):
    header = None
    seq = ""
    fasta_dict = {}

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    fasta_dict[header] = seq
                    seq = ""
                header = line[1:].split(" ")[0]
            else:
                line = line.upper()
                if not _validate_seq(line):
                    raise FastaError("Non-ACGTN charcter in \"{0}\""
                                     .format(filename))
                seq += line

        if header and len(seq):
            fasta_dict[header] = seq

    return fasta_dict


def write_fasta_dict(fasta_dict, filename):
    with open(filename, "w") as f:
        for header, seq in _iter_dict(fasta_dict):
            f.write(">{0}\n".format(header))

            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + "\n")


def reverse_complement(string):
    res = "".join(map(_comp_sym, string[::-1]))
    return res


def _validate_seq(sequence):
    for c in sequence:
        if c not in "ACGTN":
            return False
    return True


def _iter_dict(d):
    iter_d = None
    try:
        iter_d = d.iteritems()
    except AttributeError:
        iter_d = d.items()
    return iter_d


COMPL = {"A" : "T", "T" : "A", "G" : "C", "C" : "G", "N" : "N"}
def _comp_sym(char):
    return COMPL[char]
