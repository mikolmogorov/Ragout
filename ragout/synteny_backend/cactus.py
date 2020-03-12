#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module runs progressiveCactus
"""

from __future__ import absolute_import
from __future__ import division
import os
import shutil
import subprocess
import logging

from .synteny_backend import SyntenyBackend, BackendException
from .hal import HalBackend

CACTUS_EXEC = "bin/runProgressiveCactus.sh"
CACTUS_WORKDIR = "cactus-workdir"

try:
    CACTUS_INSTALL = os.environ["CACTUS_INSTALL"]
except Exception:
    CACTUS_INSTALL = ""
logger = logging.getLogger()

class CactusBackend(SyntenyBackend):
    def __init__(self):
        SyntenyBackend.__init__(self)


    def run_backend(self, recipe, output_dir, overwrite):
        logger.warning("cactus backend support is deprecated and will be "
                       "removed in future releases. Use hal instead.")
        return _make_permutations(recipe, output_dir, overwrite, self.threads)


if os.path.isfile(os.path.join(CACTUS_INSTALL, CACTUS_EXEC)):
    SyntenyBackend.register_backend("cactus", CactusBackend())


def _make_permutations(recipe, output_dir, overwrite, threads):
    """
    Runs Cactus, then outputs preprocessesed results into output_dir
    """
    work_dir = os.path.join(output_dir, CACTUS_WORKDIR)

    if overwrite and os.path.isdir(work_dir):
        shutil.rmtree(work_dir)

    hal_backend = HalBackend()
    if os.path.isdir(work_dir):
        #using existing results
        logger.warning("Using existing Cactus results from previous run")
        logger.warning("Use --overwrite to force alignment")
        hal_file = os.path.join(work_dir, "alignment.hal")
        if not os.path.isfile(hal_file):
            raise BackendException("Exitsing results are incompatible "
                                   "with input recipe")
        recipe["hal"] = hal_file

    else:
        #running cactus
        for genome, params in recipe["genomes"].items():
            if "fasta" not in params:
                raise BackendException("FASTA file for {0} is not "
                                       "specified".format(genome))

        os.mkdir(work_dir)
        config_path = _make_cactus_config(recipe, work_dir)
        hal_file = _run_cactus(config_path, work_dir, threads)
        recipe["hal"] = hal_file

    #now run another backend
    return hal_backend.run_backend(recipe, output_dir, overwrite)


def _make_cactus_config(recipe, directory):
    """
    Creates cactus "seq" file
    """
    CONF_NAME = "cactus.cfg"
    conf_file = open(os.path.join(directory, CONF_NAME), "w")
    conf_file.write(recipe["tree"] + "\n")

    for genome, params in recipe["genomes"].items():
        if genome != recipe["target"]:
            conf_file.write("*")
        conf_file.write("{0} {1}\n".format(genome, os.path.abspath(params["fasta"])))

    return conf_file.name


def _run_cactus(config_path, out_dir, threads):
    """
    Runs Progressive Cactus
    """
    CACTUS_OUT = "alignment.hal"

    logger.info("Running progressiveCactus...")
    work_dir = os.path.abspath(out_dir)
    out_hal = os.path.join(work_dir, CACTUS_OUT)
    config_file = os.path.abspath(config_path)

    threads_param = "--maxThreads=" + str(threads)

    os.chdir(CACTUS_INSTALL)
    devnull = open(os.devnull, "w")
    cmdline = [CACTUS_EXEC, config_file, work_dir, out_hal, threads_param]
    try:
        subprocess.check_call(cmdline, stdout=devnull)
    except subprocess.CalledProcessError:
        raise BackendException()

    return out_hal
