#This module runs Sibelia
##############################################################

import os
import sys
import shutil
import subprocess
import copy
import logging

from ragout.shared import utils
from .synteny_backend import SyntenyBackend, BackendException

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

    def run_backend(self, recipe, output_dir, overwrite):
        files = {}
        work_dir = os.path.join(output_dir, SIBELIA_WORKDIR)
        if overwrite and os.path.isdir(work_dir):
            shutil.rmtree(work_dir)

        if os.path.isdir(work_dir):
            #using existing results
            logger.warning("Using existing Sibelia results from previous run")
            logger.warning("Use --overwrite to force alignment")
            for block_size in recipe["blocks"]:
                block_dir = os.path.join(work_dir, str(block_size))
                perm_file = os.path.join(block_dir, "genomes_permutations.txt")
                if not os.path.isfile(perm_file):
                    raise BackendException("Exitsing results are incompatible "
                                           "with input recipe")
                files[block_size] = os.path.abspath(perm_file)

        else:
            for genome, params in recipe["genomes"].items():
                if "fasta" not in params:
                    raise BackendException("FASTA file for '{0}' is not "
                                           "specified".format(genome))

            os.mkdir(work_dir)
            chr2genome = _get_chr2genome(recipe)
            for block_size in recipe["blocks"]:
                block_dir = os.path.join(work_dir, str(block_size))
                if not os.path.isdir(block_dir):
                    os.mkdir(block_dir)

                all_fasta = [p["fasta"] for p in recipe["genomes"].values()]
                perm_file = _run_sibelia(all_fasta, block_size, block_dir)
                _postprocess(chr2genome, perm_file)
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

def _get_chr2genome(recipe):
    chr2genome = {}
    #for gen_name, fasta in genomes.items():
    for gen_name, gen_params in recipe["genomes"].items():
        with open(gen_params["fasta"], "r") as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    chr_name = line.strip()[1:].split(" ")[0]
                    chr2genome[chr_name] = gen_name
    return chr2genome


def _postprocess(chr2genome, file):
    new_file = file + "_new"
    with open(file, "r") as fin, open(new_file, "w") as fout:
        for line in fin:
            line = line.strip()
            if line.startswith(">"):
                chr_name = line[1:]
                fout.write(">{0}.{1}\n".format(chr2genome[chr_name], chr_name))
            else:
                fout.write(line + "\n")
    os.remove(file)
    os.rename(new_file, file)


def _run_sibelia(fasta_files, block_size, out_dir):

    logger.info("Running Sibelia with block size " + str(block_size))
    if not utils.which(SIBELIA_EXEC):
        raise BackendException("Sibelia is not installed")

    devnull = open(os.devnull, "w")
    cmdline = [SIBELIA_EXEC, "-s", "fine", "-m", str(block_size), "-o", out_dir]
    cmdline.extend(fasta_files)
    subprocess.check_call(cmdline, stdout=devnull)

    os.remove(os.path.join(out_dir, "d3_blocks_diagram.html"))
    shutil.rmtree(os.path.join(out_dir, "circos"))

    return os.path.join(out_dir, "genomes_permutations.txt")
