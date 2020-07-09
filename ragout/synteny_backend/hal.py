#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module works with MAF input and converts it into synteny blocks
"""

from __future__ import absolute_import
from __future__ import division
import os
import logging
import shutil
import subprocess

from .synteny_backend import SyntenyBackend, BackendException
import ragout.maf2synteny.maf2synteny as m2s
from ragout.shared import config
from ragout.shared import utils

logger = logging.getLogger()

HAL_WORKDIR = "hal-workdir"
HAL2MAF = "hal2mafMP.py"
HAL2FASTA = "hal2fasta"
HAL_STATS = "halStats"
TARGET_FASTA = "target.fasta"


class HalBackend(SyntenyBackend):
    def __init__(self):
        SyntenyBackend.__init__(self)
        self.target_fasta = None

    def infer_block_scale(self, recipe):
        hal = recipe.get("hal")
        if not hal or not os.path.exists(hal):
            raise BackendException("Could not open HAL file "
                                   "or it is not specified")
        stats = subprocess.check_output([HAL_STATS, hal])
        size = 0
        for line in stats.splitlines():
            line = line.decode()
            tokens = line.split(",")
            if tokens[0] == recipe["target"]:
                size = int(tokens[2])

        if size < config.vals["big_genome_threshold"]:
            return "small"
        else:
            return "large"

    def run_backend(self, recipe, output_dir, overwrite):
        workdir = os.path.join(output_dir, HAL_WORKDIR)
        if overwrite and os.path.isdir(workdir):
            shutil.rmtree(workdir)

        if "hal" not in recipe or not os.path.exists(recipe["hal"]):
            raise BackendException("Could not open HAL file "
                                   "or it is not specified")

        files = {}
        #using existing results
        if os.path.isdir(workdir):
            logger.warning("Using synteny blocks from previous run")
            logger.warning("Use --overwrite to force alignment")

            all_good = True
            for block_size in self.blocks:
                block_dir = os.path.join(workdir, str(block_size))
                coords_file = os.path.join(block_dir, "blocks_coords.txt")
                if not os.path.isfile(coords_file):
                    all_good = False
                    break
                files[block_size] = os.path.abspath(coords_file)

            target_fasta = os.path.join(workdir, TARGET_FASTA)
            if not os.path.isfile(target_fasta):
                all_good = False
            else:
                self.target_fasta = target_fasta

            if not all_good:
                raise BackendException("Exitsing results are incompatible "
                                           "with current run")

        else:
            os.mkdir(workdir)

            logger.info("Extracting FASTA from HAL")
            target_fasta = os.path.join(workdir, TARGET_FASTA)
            cmdline = [HAL2FASTA, recipe["hal"], recipe["target"],
                       "--inMemory"]
            subprocess.check_call(cmdline, stdout=open(target_fasta, "w"))
            self.target_fasta = target_fasta

            logger.info("Converting HAL to MAF")
            out_maf = os.path.join(workdir, "alignment.maf")
            ref_genome = recipe["target"]   #Tricky notation, huh?
            export_genomes = ",".join(recipe["genomes"])

            cmdline = [HAL2MAF, recipe["hal"], out_maf, "--noAncestors",
                        "--numProc", str(self.threads),  "--refGenome",
                        ref_genome, "--targetGenomes", export_genomes]
            logger.debug(" ".join(cmdline))
            subprocess.check_call(cmdline, stdout=open(os.devnull, "w"))

            logger.info("Extracting synteny blocks from MAF")
            if not m2s.make_synteny(out_maf, workdir, self.blocks):
                raise BackendException("Something went wrong with maf2synteny")

            for block_size in self.blocks:
                block_dir = os.path.join(workdir, str(block_size))
                coords_file = os.path.join(block_dir, "blocks_coords.txt")
                files[block_size] = os.path.abspath(coords_file)
                if not os.path.exists(coords_file):
                    raise BackendException("Something bad happened!")

        return files


if utils.which(HAL2MAF) and utils.which(HAL2FASTA) and utils.which(HAL_STATS):
    SyntenyBackend.register_backend("hal", HalBackend())
