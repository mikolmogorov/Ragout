#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module works with MAF input and converts it into synteny blocks
"""

import os
import logging
import shutil
import multiprocessing
import subprocess

from .synteny_backend import SyntenyBackend, BackendException
import ragout.maf2synteny.maf2synteny as m2s
from ragout.shared import config
from ragout.shared import utils

logger = logging.getLogger()
HAL_WORKDIR = "hal-workdir"
HAL2MAF = "hal2mafMP.py"


class HalBackend(SyntenyBackend):
    def __init__(self):
        SyntenyBackend.__init__(self)

    def run_backend(self, recipe, output_dir, overwrite):
        workdir = os.path.join(output_dir, HAL_WORKDIR)
        if overwrite and os.path.isdir(workdir):
            shutil.rmtree(workdir)

        if "hal" not in recipe or not os.path.exists(recipe["hal"]):
            raise BackendException("Could not open HAL file "
                                   "or it is not specified")

        files = {}
        if os.path.isdir(workdir):
            #using existing results
            logger.warning("Using synteny blocks from previous run")
            logger.warning("Use --overwrite to force alignment")
            for block_size in recipe["blocks"]:
                block_dir = os.path.join(workdir, str(block_size))
                coords_file = os.path.join(block_dir, "blocks_coords.txt")
                if not os.path.isfile(coords_file):
                    raise BackendException("Exitsing results are incompatible "
                                           "with input recipe")
                files[block_size] = os.path.abspath(coords_file)

        else:
            os.mkdir(workdir)
            ###Running hal2maf
            logger.info("Converting HAL to MAF")
            num_proc = min(config.vals["cactus_max_threads"],
                           multiprocessing.cpu_count())
            out_maf = os.path.join(workdir, "alignment.maf")
            ref_genome = recipe["target"]   #Tricky notation, huh?
            export_genomes = ",".join(recipe["genomes"])

            cmdline = [HAL2MAF, recipe["hal"], out_maf, "--noAncestors",
                        "--numProc", str(num_proc),  "--refGenome", ref_genome,
                        "--targetGenomes", export_genomes, "--inMemory"]
            logger.debug(" ".join(cmdline))
            subprocess.check_call(cmdline, stdout=open(os.devnull, "w"))

            ###Running maf2synteny
            logger.info("Extracting synteny blocks from MAF")
            if not m2s.make_synteny(recipe["maf"], workdir, recipe["blocks"]):
                raise BackendException("Something went wrong with maf2synteny")

            for block_size in recipe["blocks"]:
                block_dir = os.path.join(workdir, str(block_size))
                coords_file = os.path.join(block_dir, "blocks_coords.txt")
                files[block_size] = os.path.abspath(coords_file)
                if not os.path.exists(coords_file):
                    raise BackendException("Something bad happened!")

        return files


if utils.which(HAL2MAF):
    SyntenyBackend.register_backend("hal", HalBackend())
