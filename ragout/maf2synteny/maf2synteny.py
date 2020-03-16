#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module executes maf2synteny native binary
which recovers synteny blocks from multiple alignment
"""

from __future__ import absolute_import
from __future__ import division
import logging
import subprocess
import os

from ragout.shared.utils import which
from ragout.shared import config
from ragout.six.moves import map

logger = logging.getLogger()

M2S_EXEC = "ragout-maf2synteny"


def check_binary():
    """
    Checks if native binary is available and runnable
    """
    binary = which(M2S_EXEC)
    if not binary:
        logger.error("\"%s\" native module not found", M2S_EXEC)
        return False

    try:
        devnull = open(os.devnull, "w")
        subprocess.check_call([M2S_EXEC, "-h"], stderr=devnull)
    except subprocess.CalledProcessError as e:
        logger.error("Some error inside native module: %s", str(e))
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
    cmdline = [M2S_EXEC, maf_file, "-o", out_dir, "-s", params_file,
               "-b", ",".join(map(str, min_blocks_list))]

    logger.info("Running maf2synteny module")
    proc = subprocess.Popen(cmdline)
    #for line in iter(proc.stderr.readline, ""):
    #    logger.debug(line.strip())
    ret_code = proc.wait()
    if ret_code:
        logger.error("Non-zero return code: %d", ret_code)
        return False

    os.remove(params_file)

    return True

def _make_params_file(params, out_file):
    assert len(params)

    with open(out_file, "w") as f:
        for k, d in params:
            f.write("{0} {1}\n".format(k, d))
