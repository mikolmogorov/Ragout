#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module executes overlap native binary
which reconstructs overlap graph from contigs
"""

from __future__ import absolute_import
from __future__ import division
import os
import logging
import subprocess

from ragout.shared import config
from ragout.shared.utils import which

logger = logging.getLogger()

OVERLAP_EXEC = "ragout-overlap"

class OverlapException(Exception):
    pass

def check_binary():
    """
    Checks if the native binary is available and runnable
    """
    binary = which(OVERLAP_EXEC)
    if not binary:
        logger.error("\"%s\" native module not found", OVERLAP_EXEC)
        return False

    try:
        devnull = open(os.devnull, "w")
        subprocess.check_call([OVERLAP_EXEC, "--help"], stderr=devnull)
    except subprocess.CalledProcessError as e:
        logger.error("Some error inside native module: %s", str(e))
        return False

    return True


def make_overlap_graph(contigs_file, dot_file):
    """
    Builds assembly graph and outputs it in "dot" format
    """
    cmdline = [OVERLAP_EXEC, contigs_file, dot_file,
               str(config.vals["overlap"]["min_overlap"]),
               str(config.vals["overlap"]["max_overlap"])]
    if config.vals["overlap"]["detect_kmer"]:
        cmdline.append("--detect-kmer")

    logger.info("Building assembly graph")
    proc = subprocess.Popen(cmdline)
    #for line in iter(proc.stderr.readline, ""):
    #    logger.debug(line.strip())
    ret_code = proc.wait()
    if ret_code:
        raise OverlapException("Error building overlap graph: Non-zero return "
                               "code when calling {0} "
                               "module: {1}".format(OVERLAP_EXEC, ret_code))
