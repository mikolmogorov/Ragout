#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module executes maf2synteny native binary
which recovers synteny blocks from multiple alignment
"""

import logging
import subprocess
import os

from ragout.shared.utils import which
from ragout.shared import config

logger = logging.getLogger()

M2S_EXEC = "ragout-maf2synteny"


def check_binary():
    """
    Checks if native binary is available and runnable
    """
    binary = which(M2S_EXEC)
    if not binary:
        logger.error("\"{0}\" native module not found".format(M2S_EXEC))
        return False

    try:
        devnull = open(os.devnull, "w")
        subprocess.check_call([M2S_EXEC, "--help"], stderr=devnull)
    except subprocess.CalledProcessError as e:
        logger.error("Some error inside native {0} module: {1}"
                     .format(M2S_EXEC, e))
        return False

    return True


def make_synteny(maf_file, out_dir, min_blocks_list):
    """
    Builds synteny blocks from MAF file
    """
    if not check_binary():
        return False

    params_file = os.path.join(out_dir, "simpl_params.txt")
    _make_params_file(config.vals["maf2synteny"], params_file)
    cmdline = [M2S_EXEC, maf_file, out_dir, params_file]
    cmdline.extend(list(map(str, min_blocks_list)))

    logger.info("Running maf2synteny module")
    proc = subprocess.Popen(cmdline, stderr=subprocess.PIPE)
    for line in iter(proc.stderr.readline, ""):
        logger.debug(line.strip())
    ret_code = proc.wait()
    if ret_code:
        logger.error("Non-zero return code when calling {0} module: {1}"
                     .format(M2S_EXEC, ret_code))
        return False

    os.remove(params_file)

    return True

def _make_params_file(params, out_file):
    assert len(params)

    with open(out_file, "w") as f:
        for k, d in params:
            f.write("{0} {1}\n".format(k, d))
