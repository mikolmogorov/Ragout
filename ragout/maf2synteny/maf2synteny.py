"""
This module executes maf2synteny native binary
which recovers synteny blocks from multiple alignment
"""

import logging
import subprocess

from ragout.shared.utils import which

M2S_EXEC = "ragout-maf2synteny"


def make_synteny(maf_file, out_dir, min_blocks_list):
    """
    Builds synteny blocks from MAF file
    """
    if not which(M2S_EXEC):
        logger.error("\"{0}\" native module not found".format(M2S_EXEC))
        return False

    cmdline = [M2S_EXEC, maf_file, out_dir]
    cmdline.extend(list(map(str, min_blocks_list)))
    try:
        subprocess.check_call(cmdline)
    except subprocess.CalledProcessError as e:
        logger.error("Some error inside native {0} module".format(M2S_EXEC))
        return False

    return True
