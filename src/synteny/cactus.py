###############################################################
#This module runs progressiveCactus
###############################################################

import os
import sys
import shutil
import subprocess
import logging

from .. import utils
from synteny_backend import SyntenyBackend
import maf2synteny.maf2synteny as m2s

logger = logging.getLogger()


#PUBLIC:
################################################################

class CactusBackend(SyntenyBackend):
    def __init__(self):
        SyntenyBackend.__init__(self)


    def run_backend(self, config, output_dir):
        return make_permutations(config.references, config.targets,
                                 config.tree, config.blocks, output_dir)


if True:
    logger.debug("progressiveCactus is installed")
    SyntenyBackend.register_backend("cactus", CactusBackend())
else:
    logger.debug("progressiveCactus is not installed")


#PRIVATE:
#################################################################


#Runs Cactus, then outputs preprocessesed results into output_dir
def make_permutations(references, targets, tree, block_sizes, output_dir):
    config_path = make_cactus_config(references, targets, tree, output_dir)
    ref_genome = targets.keys()[0]
    maf_file = run_cactus(config_path, ref_genome, output_dir)

    files = {}
    for block_size in block_sizes:
        block_dir = os.path.join(output_dir, str(block_size))
        if not os.path.isdir(block_dir):
            os.mkdir(block_dir)
        m2s.get_synteny(maf_file, block_dir, block_size)
        perm_file = os.path.join(block_dir, "genomes_permutations.txt")

        files[block_size] = os.path.abspath(perm_file)

    return files


def make_cactus_config(references, targets, tree_string, directory):
    CONF_NAME = "cactus.cfg"
    file = open(os.path.join(directory, CONF_NAME), "w")
    file.write(tree_string + "\n")

    genomes = dict(references.items() + targets.items())
    for seq_id, seq_path in genomes.iteritems():
        file.write("{0} {1}\n".format(seq_id, os.path.abspath(seq_path)))

    return file.name


def run_cactus(config_path, ref_genome, out_dir):
    CACTUS_DIR = "/home/volrath/Bioinf/Tools/progressiveCactus/"
    CACTUS_EXEC = "bin/runProgressiveCactus.sh"
    CACTUS_OUT = "alignment.hal"
    HAL2MAF = "submodules/hal/bin/hal2maf"
    MAF_OUT = "cactus.maf"

    logger.info("Running progressiveCactus...")
    work_dir = os.path.abspath(os.path.join(out_dir, "cactus-workdir"))
    out_hal = os.path.abspath(os.path.join(work_dir, CACTUS_OUT))
    out_maf = os.path.abspath(os.path.join(out_dir, MAF_OUT))
    config_file = os.path.abspath(config_path)
    prev_dir = os.getcwd()
    if not os.path.exists(CACTUS_DIR):
        raise Exception("progressiveCactus is not installed")

    os.chdir(CACTUS_DIR)
    #devnull = open(os.devnull, "w")
    cmdline = [CACTUS_EXEC, config_file, work_dir, out_hal]
    subprocess.check_call(cmdline)

    #convert to maf
    logger.info("Converting HAL to MAF...")
    cmdline = [HAL2MAF, out_hal, out_maf, "--noAncestors", "--refGenome", ref_genome]
    subprocess.check_call(cmdline)

    os.chdir(prev_dir)
    return out_maf


