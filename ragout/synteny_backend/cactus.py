###############################################################
#This module runs progressiveCactus
###############################################################

import os
import sys
import shutil
import subprocess
import multiprocessing
import logging

from .synteny_backend import SyntenyBackend
import ragout.maf2synteny.maf2synteny as m2s

CACTUS_EXEC = "bin/runProgressiveCactus.sh"
CACTUS_WORKDIR = "cactus-workdir"
try:
    CACTUS_INSTALL = os.environ["CACTUS_INSTALL"]
except:
    CACTUS_INSTALL = ""
logger = logging.getLogger()

#PUBLIC:
################################################################

class CactusBackend(SyntenyBackend):
    def __init__(self):
        SyntenyBackend.__init__(self)


    def run_backend(self, config, output_dir, overwrite):
        return _make_permutations(config.references, config.targets, config.tree,
                                  config.blocks, output_dir, overwrite)


if os.path.isfile(os.path.join(CACTUS_INSTALL, CACTUS_EXEC)):
    logger.debug("progressiveCactus is installed")
    SyntenyBackend.register_backend("cactus", CactusBackend())
else:
    logger.debug("progressiveCactus is not installed")


#PRIVATE:
#################################################################


#Runs Cactus, then outputs preprocessesed results into output_dir
def _make_permutations(references, targets, tree, block_sizes,
                       output_dir, overwrite):
    work_dir = os.path.join(output_dir, CACTUS_WORKDIR)
    files = {}

    if overwrite and os.path.isdir(work_dir):
        shutil.rmtree(work_dir)

    if os.path.isdir(work_dir):
        #using existing results
        logger.warning("Using existing Cactus results from previous run")
        logger.warning("Use --overwrite to force alignment")
        for block_size in block_sizes:
            block_dir = os.path.join(work_dir, str(block_size))
            perm_file = os.path.join(block_dir, "genomes_permutations.txt")
            if not os.path.isfile(perm_file):
                logger.error("Exitsing results are incompatible with input config")
                raise Exception("Cannot reuse results from previous run")
            files[block_size] = os.path.abspath(perm_file)

    else:
        #running cactus
        os.mkdir(work_dir)
        config_path = _make_cactus_config(references, targets, tree, work_dir)
        ref_genome = targets.keys()[0]
        maf_file = _run_cactus(config_path, ref_genome, work_dir)

        for block_size in block_sizes:
            block_dir = os.path.join(work_dir, str(block_size))
            if not os.path.isdir(block_dir):
                os.mkdir(block_dir)
            if not m2s.make_synteny(maf_file, block_dir, block_size):
                raise Exception("Something went wrong with maf2synteny")

            perm_file = os.path.join(block_dir, "genomes_permutations.txt")
            files[block_size] = os.path.abspath(perm_file)

    return files


def _make_cactus_config(references, targets, tree_string, directory):
    CONF_NAME = "cactus.cfg"
    file = open(os.path.join(directory, CONF_NAME), "w")
    file.write(tree_string + "\n")

    #genomes = dict(references.items() + targets.items())
    for seq_id, seq_path in references.iteritems():
        file.write("*{0} {1}\n".format(seq_id, os.path.abspath(seq_path)))
    for seq_id, seq_path in targets.iteritems():
        file.write("{0} {1}\n".format(seq_id, os.path.abspath(seq_path)))

    return file.name


def _run_cactus(config_path, ref_genome, out_dir):
    CACTUS_OUT = "alignment.hal"
    HAL2MAF = "submodules/hal/bin/hal2maf"
    MAF_OUT = "cactus.maf"
    MAX_THREADS = 10

    logger.info("Running progressiveCactus...")
    work_dir = os.path.abspath(out_dir)
    out_hal = os.path.join(work_dir, CACTUS_OUT)
    out_maf = os.path.join(work_dir, MAF_OUT)
    config_file = os.path.abspath(config_path)
    prev_dir = os.getcwd()
    #if not os.path.exists(CACTUS_DIR):
    #    raise Exception("progressiveCactus is not installed")

    num_proc = min(MAX_THREADS, multiprocessing.cpu_count())
    threads_param = "--maxThreads=" + str(num_proc)

    os.chdir(CACTUS_INSTALL)
    devnull = open(os.devnull, "w")
    cmdline = [CACTUS_EXEC, config_file, work_dir, out_hal, threads_param]
    subprocess.check_call(cmdline, stdout=devnull)

    #convert to maf
    logger.info("Converting HAL to MAF...")
    cmdline = [HAL2MAF, out_hal, out_maf, "--noAncestors",
               "--refGenome", ref_genome]
    subprocess.check_call(cmdline, stdout=devnull)

    os.chdir(prev_dir)
    return out_maf
