#This module runs Sibelia
##############################################################

import os
import sys
import shutil
import subprocess
import copy
import logging

from ragout.shared import utils
from .synteny_backend import SyntenyBackend

logger = logging.getLogger()

SIBELIA_EXEC = "Sibelia"
SIBELIA_WORKDIR = "sibelia-workdir"
try:
    SIBELIA_INSTALL = os.environ["SIBELIA_INSTALL"]
    os.environ["PATH"] += os.pathsep + SIBELIA_INSTALL
except:
    pass


#PUBLIC:
################################################################

class SibeliaBackend(SyntenyBackend):
    def __init__(self):
        SyntenyBackend.__init__(self)

    def run_backend(self, config, output_dir, overwrite):
        files = {}
        work_dir = os.path.join(output_dir, SIBELIA_WORKDIR)
        if overwrite and os.path.isdir(work_dir):
            shutil.rmtree(work_dir)

        if os.path.isdir(work_dir):
            #using existing results
            logger.warning("Using existing Sibelia results from previous run")
            logger.warning("Use --overwrite to force alignment")
            for block_size in config.blocks:
                block_dir = os.path.join(work_dir, str(block_size))
                perm_file = os.path.join(block_dir, "genomes_permutations.txt")
                if not os.path.isfile(perm_file):
                    logger.error("Exitsing results are incompatible with input config")
                    raise Exception
                files[block_size] = os.path.abspath(perm_file)

        else:
            os.mkdir(work_dir)
            genomes = dict(config.references.items() + config.targets.items())
            for block_size in config.blocks:
                block_dir = os.path.join(work_dir, str(block_size))
                if not os.path.isdir(block_dir):
                    os.mkdir(block_dir)

                perm_file = _run_sibelia(genomes.values(), block_size, block_dir)
                files[block_size] = perm_file

        return files


def _check_installation():
    return bool(utils.which(SIBELIA_EXEC))

if _check_installation():
    logger.debug("Sibelia is installed")
    SyntenyBackend.register_backend("sibelia", SibeliaBackend())
else:
    logger.debug("Sibelia is not installed")


#PRIVATE:
#################################################################

def _run_sibelia(fasta_files, block_size, out_dir):

    logger.info("Running Sibelia with block size " + str(block_size))
    if not utils.which(SIBELIA_EXEC):
        raise Exception("Sibelia is not installed")

    devnull = open(os.devnull, "w")
    cmdline = [SIBELIA_EXEC, "-s", "fine", "-m", str(block_size), "-o", out_dir]
    cmdline.extend(fasta_files)
    subprocess.check_call(cmdline, stdout=devnull)

    #os.remove(os.path.join(out_dir, "coverage_report.txt"))
    #os.remove(os.path.join(out_dir, "blocks_coords.txt"))
    os.remove(os.path.join(out_dir, "d3_blocks_diagram.html"))
    shutil.rmtree(os.path.join(out_dir, "circos"))

    return os.path.join(out_dir, "genomes_permutations.txt")
