#This module runs Sibelia
##############################################################

import os
import sys
import shutil
import subprocess
import copy
import logging

from .. import utils
from synteny_backend import SyntenyBackend

logger = logging.getLogger()

LIB_DIR = "lib"
SIBELIA_DIR = os.path.join(LIB_DIR, "Sibelia")
#running_dir = os.path.dirname(os.path.realpath(__file__))
running_dir = os.getcwd()
os.environ["PATH"] += os.pathsep + os.path.join(running_dir, SIBELIA_DIR)
SIBELIA_EXEC = "Sibelia"


#PUBLIC:
################################################################

class SibeliaBackend(SyntenyBackend):
    def __init__(self):
        SyntenyBackend.__init__(self)


    def run_backend(self, config, output_dir):
        files = {}
        genomes = dict(config.references.items() + config.targets.items())

        for block_size in config.blocks:
            block_dir = os.path.join(output_dir, str(block_size))
            if not os.path.isdir(block_dir):
                os.mkdir(block_dir)

            perm_file = run_sibelia(genomes.values(), block_size, block_dir)
            files[block_size] = perm_file

        return files


def check_installation():
    return bool(utils.which(SIBELIA_EXEC))

if check_installation():
    logger.debug("Sibelia is installed")
    SyntenyBackend.register_backend("sibelia", SibeliaBackend())
else:
    logger.debug("Sibelia is not installed")


#PRIVATE:
#################################################################

def run_sibelia(fasta_files, block_size, out_dir):

    logger.info("Running Sibelia with block size " + str(block_size))
    if not utils.which(SIBELIA_EXEC):
        raise Exception("Sibelia is not installed")

    devnull = open(os.devnull, "w")
    cmdline = [SIBELIA_EXEC, "-s", "loose", "-m", str(block_size), "-o", out_dir]
    cmdline.extend(fasta_files)
    subprocess.check_call(cmdline, stdout=devnull)

    os.remove(os.path.join(out_dir, "coverage_report.txt"))
    os.remove(os.path.join(out_dir, "d3_blocks_diagram.html"))
    os.remove(os.path.join(out_dir, "blocks_coords.txt"))
    shutil.rmtree(os.path.join(out_dir, "circos"))

    return os.path.join(out_dir, "genomes_permutations.txt")
