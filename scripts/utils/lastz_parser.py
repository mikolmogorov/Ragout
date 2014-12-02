
#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
Functions for lastz input handling
"""

import subprocess
import os

from .common import AlignmentInfo

def parse_lastz_maf(filename):
    alignments = []
    with open(filename, "r") as f:
        while True:
            line = f.readline()
            if not line:
                break

            if not line.startswith("a"):
                state = 1
                continue

            #read two next lines
            ref_line = f.readline()
            qry_line = f.readline()
            assert ref_line.startswith("s") and qry_line.startswith("s")

            ref_vals = ref_line.strip().split()
            ref_name, ref_start, ref_size, ref_strand, ref_len = ref_vals[1:6]
            ref_len, ref_size = int(ref_len), int(ref_size)
            if ref_strand == "+":
                ref_start, ref_end = int(ref_start), int(ref_start) + ref_size
            else:
                ref_start = ref_len - 1 - int(ref_start)
                ref_end = ref_start - ref_size

            qry_vals = qry_line.strip().split()
            qry_name, qry_start, qry_size, qry_strand, qry_len = qry_vals[1:6]
            qry_len, qry_size = int(ref_len), int(ref_size)
            if qry_strand == "+":
                qry_start, qry_end = int(qry_start), int(qry_start) + qry_size
            else:
                qry_start = qry_len - 1 - int(qry_start)
                qry_end = qry_start - qry_size

            alignments.append(AlignmentInfo(ref_start, ref_end, qry_start,
                                            qry_end, ref_len, qry_len,
                                            ref_name, qry_name))

    return alignments


LASTZ_BIN = "lastz"
def run_lastz(reference, target, out_file):
    print("Running lastz")
    reference = os.path.abspath(reference)
    target = os.path.abspath(target)
    out_file = os.path.abspath(out_file)
    cmdline = [LASTZ_BIN, reference + "[multiple,nameparse=darkspace]",
               target + "[nameparse=darkspace]","--notransition",
               "--step=20", "--chain", "--gapped", "--gfextend",
               "--ambiguous=n", "--format=maf", "--output=" + out_file]
    devnull = os.devnull
    subprocess.check_call(cmdline, stderr=open(devnull, "w"))
