#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module provides some basic FASTA I/O
"""

from string import maketrans

class FastaError(Exception):
    pass

def read_fasta_dict(filename):
    """
    Reads fasta file into dictionary. Also preforms some validation
    """
    header = None
    seq = []
    fasta_dict = {}

    try:
        with open(filename, "r") as f:
            for lineno, line in enumerate(f):
                line = line.strip()
                if line.startswith(">"):
                    if header:
                        fasta_dict[header] = "".join(seq)
                        seq = []
                    header = line[1:].split(" ")[0]
                else:
                    if not _validate_seq(line):
                        raise FastaError("Invalid char in \"{0}\" at line {1}"
                                         .format(filename, lineno))
                    seq.append(line)

            if header and len(seq):
                fasta_dict[header] = "".join(seq)

    except IOError as e:
        raise FastaError(e)

    return fasta_dict


def write_fasta_dict(fasta_dict, filename):
    """
    Writes dictionary with fasta to file
    """
    with open(filename, "w") as f:
        for header, seq in _iter_dict(fasta_dict):
            f.write(">{0}\n".format(header))

            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + "\n")


COMPL = maketrans("ATGCURYKMSWBVDHNXatgcurykmswbvdhnx",
                  "TACGAYRMKSWVBHDNXtacgayrmkswvbhdnx")
def reverse_complement(string):
    return string[::-1].translate(COMPL)


def _validate_seq(sequence):
    VALID_CHARS = "ACGTURYKMSWBDHVNXatgcurykmswbvdhnx"
    if len(sequence.translate(None, VALID_CHARS)):
        return False
    return True


def _iter_dict(d):
    """
    Helper fundtion to iterate through dictionary in both Python 2 and 3
    """
    iter_d = None
    try:
        iter_d = d.iteritems()
    except AttributeError:
        iter_d = d.items()
    return iter_d
