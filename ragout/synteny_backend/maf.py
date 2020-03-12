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

from .synteny_backend import SyntenyBackend, BackendException
import ragout.maf2synteny.maf2synteny as m2s

logger = logging.getLogger()
MAF_WORKDIR = "maf-workdir"

class MafBackend(SyntenyBackend):
    def __init__(self):
        SyntenyBackend.__init__(self)

    def run_backend(self, recipe, output_dir, overwrite):
        logger.warning("Maf support is deprecated and will be removed "
                       "in future releases. Use hal istead.")
        workdir = os.path.join(output_dir, MAF_WORKDIR)
        if overwrite and os.path.isdir(workdir):
            shutil.rmtree(workdir)

        if "maf" not in recipe or not os.path.exists(recipe["maf"]):
            raise BackendException("Could not open MAF file "
                                   "or it is not specified")

        files = {}
        if os.path.isdir(workdir):
            #using existing results
            logger.warning("Using synteny blocks from previous run")
            logger.warning("Use --overwrite to force alignment")
            for block_size in self.blocks:
                block_dir = os.path.join(workdir, str(block_size))
                coords_file = os.path.join(block_dir, "blocks_coords.txt")
                if not os.path.isfile(coords_file):
                    raise BackendException("Exitsing results are incompatible "
                                           "with input recipe")
                files[block_size] = os.path.abspath(coords_file)

        else:
            os.mkdir(workdir)
            logger.info("Converting MAF to synteny")
            if not m2s.make_synteny(recipe["maf"], workdir, self.blocks):
                raise BackendException("Something went wrong with maf2synteny")

            for block_size in self.blocks:
                block_dir = os.path.join(workdir, str(block_size))
                coords_file = os.path.join(block_dir, "blocks_coords.txt")
                files[block_size] = os.path.abspath(coords_file)
                if not os.path.exists(coords_file):
                    raise BackendException("Something bad happened!")

        return files

SyntenyBackend.register_backend("maf", MafBackend())
