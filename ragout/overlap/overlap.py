#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module executes overlap native binary
which reconstructs overlap graph from contigs
"""

import logging
import subprocess

from ragout.shared import config
from ragout.shared.utils import which

logger = logging.getLogger()

OVERLAP_EXEC = "ragout-overlap"

def check_binary():
    """
    Checks if the native binary is available
    """
    binary = which(OVERLAP_EXEC)
    if not binary:
        logger.error("\"{0}\" native module not found".format(OVERLAP_EXEC))
        return False
    return True


def make_overlap_graph(contigs_file, dot_file):
    """
    Builds assembly graph and outputs it in "dot" format
    """
    logger.info("Building overlap graph...")
    if not check_binary():
        return False

    cmdline = [OVERLAP_EXEC, contigs_file, dot_file,
               str(config.vals["overlap"]["min_overlap"]),
               str(config.vals["overlap"]["max_overlap"])]
    try:
        subprocess.check_call(cmdline)
    except subprocess.CalledProcessError as e:
        logger.error("Some error inside native {0} module".format(OVERLAP_EXEC))
        return False

    return True
