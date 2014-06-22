#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module executes maf2synteny native binary
which recovers synteny blocks from multiple alignment
"""

import logging
import subprocess

from ragout.shared.utils import which

logger = logging.getLogger()

M2S_EXEC = "ragout-maf2synteny"


def check_binary():
    """
    Checks if native binary is available
    """
    binary = which(M2S_EXEC)
    if not binary:
        logger.error("\"{0}\" native module not found".format(M2S_EXEC))
        return False
    return True


def make_synteny(maf_file, out_dir, min_blocks_list):
    """
    Builds synteny blocks from MAF file
    """
    if not check_binary():
        return False

    cmdline = [M2S_EXEC, maf_file, out_dir]
    cmdline.extend(list(map(str, min_blocks_list)))
    try:
        subprocess.check_call(cmdline)
    except subprocess.CalledProcessError as e:
        logger.error("Some error inside native {0} module".format(M2S_EXEC))
        return False

    return True
