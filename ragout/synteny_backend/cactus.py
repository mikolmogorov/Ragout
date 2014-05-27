###############################################################
#This module runs progressiveCactus
###############################################################

import os
import sys
import shutil
import subprocess
import multiprocessing
import logging

from .synteny_backend import SyntenyBackend, BackendException
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


    def run_backend(self, recipe, output_dir, overwrite):
        return _make_permutations(recipe, output_dir, overwrite)


if os.path.isfile(os.path.join(CACTUS_INSTALL, CACTUS_EXEC)):
    logger.debug("progressiveCactus is installed")
    SyntenyBackend.register_backend("cactus", CactusBackend())
else:
    logger.debug("progressiveCactus is not installed")


#PRIVATE:
#################################################################


#Runs Cactus, then outputs preprocessesed results into output_dir
def _make_permutations(recipe, output_dir, overwrite):
    work_dir = os.path.join(output_dir, CACTUS_WORKDIR)
    files = {}

    if overwrite and os.path.isdir(work_dir):
        shutil.rmtree(work_dir)

    if os.path.isdir(work_dir):
        #using existing results
        logger.warning("Using existing Cactus results from previous run")
        logger.warning("Use --overwrite to force alignment")
        for block_size in recipe["blocks"]:
            block_dir = os.path.join(work_dir, str(block_size))
            perm_file = os.path.join(block_dir, "genomes_permutations.txt")
            if not os.path.isfile(perm_file):
                raise BackendException("Exitsing results are incompatible "
                                       "with input recipe")
            files[block_size] = os.path.abspath(perm_file)

    else:
        #running cactus
        for genome, params in recipe["genomes"].items():
            if "fasta" not in params:
                raise BackendException("FASTA file for {0} is not "
                                       "specified".format(genome))

        os.mkdir(work_dir)
        config_path = _make_cactus_config(recipe, work_dir)
        ref_genome = recipe["target"][0]
        maf_file = _run_cactus(config_path, ref_genome, work_dir)

        logger.info("Converting maf to synteny")
        if not m2s.make_synteny(maf_file, work_dir, recipe["blocks"]):
            raise BackendException("Something went wrong with maf2synteny")

        for block_size in recipe["blocks"]:
            block_dir = os.path.join(work_dir, str(block_size))
            perm_file = os.path.join(block_dir, "genomes_permutations.txt")
            files[block_size] = os.path.abspath(perm_file)
            if not os.path.exists(perm_file):
                raise BackendException("Something bad happened!")

    return files


def _make_cactus_config(recipe, directory):
    CONF_NAME = "cactus.cfg"
    file = open(os.path.join(directory, CONF_NAME), "w")
    file.write(recipe["tree"] + "\n")

    for genome, params in recipe["genomes"].items():
        if genome not in recipe["target"]:
            file.write("*")
        file.write("{0} {1}\n".format(genome, os.path.abspath(params["fasta"])))

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

    num_proc = min(MAX_THREADS, multiprocessing.cpu_count())
    threads_param = "--maxThreads=" + str(num_proc)

    os.chdir(CACTUS_INSTALL)
    devnull = open(os.devnull, "w")
    cmdline = [CACTUS_EXEC, config_file, work_dir, out_hal, threads_param]
    try:
        subprocess.check_call(cmdline, stdout=devnull)
    except subprocess.CalledProcessError:
        raise BackendException()

    #convert to maf
    logger.info("Converting HAL to MAF...")
    cmdline = [HAL2MAF, out_hal, out_maf, "--noAncestors",
               "--refGenome", ref_genome, "--ucscNames"]
    subprocess.check_call(cmdline, stdout=devnull)

    os.chdir(prev_dir)
    return out_maf
