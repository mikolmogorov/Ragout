#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module provides some basic FASTA I/O
"""

from __future__ import absolute_import
from __future__ import division
import logging
import sys
if sys.version_info < (3, 0):
    from string import maketrans
    _STR = lambda x: x
    _BYTES = lambda x: x
else:
    maketrans = bytes.maketrans
    _STR = bytes.decode
    _BYTES = str.encode

from ragout.six.moves import range

logger = logging.getLogger()

class FastaError(Exception):
    pass

def read_fasta_dict(filename):
    """
    Reads fasta file into dictionary. Also preforms some validation
    """
    logger.info("Reading contigs file")

    header = None
    seq = []
    fasta_dict = {}

    try:
        with open(filename, "rb") as f:
            for lineno, line in enumerate(f):
                line = line.strip()
                if line.startswith(b">"):
                    if header:
                        fasta_dict[_STR(header)] = _STR(b"".join(seq))
                        seq = []
                    header = line[1:].split(b" ")[0]
                else:
                    if not _validate_seq(line):
                        raise FastaError("Invalid char in \"{0}\" at line {1}"
                                         .format(filename, lineno))
                    seq.append(line)

            if header and len(seq):
                fasta_dict[_STR(header)] = _STR(b"".join(seq))

    except IOError as e:
        raise FastaError(e)

    return fasta_dict


def write_fasta_dict(fasta_dict, filename):
    """
    Writes dictionary with fasta to file
    """
    with open(filename, "w") as f:
        for header in sorted(fasta_dict):
            f.write(">{0}\n".format(header))

            for i in range(0, len(fasta_dict[header]), 60):
                f.write(fasta_dict[header][i:i + 60] + "\n")


COMPL = maketrans(b"ATGCURYKMSWBVDHNXatgcurykmswbvdhnx",
                  b"TACGAYRMKSWVBHDNXtacgayrmkswvbhdnx")
def reverse_complement(string):
    return _STR(_BYTES(string)[::-1].translate(COMPL))


def _validate_seq(sequence):
    VALID_CHARS = b"ACGTURYKMSWBDHVNXatgcurykmswbvdhnx"
    if len(sequence.translate(None, VALID_CHARS)):
        return False
    return True
