#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module runs Sibelia
"""

from __future__ import absolute_import
from __future__ import division
import os
import shutil
import subprocess
import logging

from ragout.shared import utils
from ragout.shared import config
from .synteny_backend import SyntenyBackend, BackendException

logger = logging.getLogger()

SIBELIA_EXEC = "Sibelia"
SIBELIA_WORKDIR = "sibelia-workdir"
SIBELIA_MAX_INPUT = 100 * 1024 * 1024

try:
    SIBELIA_INSTALL = os.environ["SIBELIA_INSTALL"]
    os.environ["PATH"] += os.pathsep + SIBELIA_INSTALL
except Exception:
    pass


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
            for block_size in self.blocks:
                block_dir = os.path.join(work_dir, str(block_size))
                coords_file = os.path.join(block_dir, "blocks_coords.txt")
                if not os.path.isfile(coords_file):
                    raise BackendException("Exitsing results are incompatible "
                                           "with input recipe")
                files[block_size] = os.path.abspath(coords_file)

        else:
            chr2genome, total_size = _get_sequence_info(recipe)
            if total_size > SIBELIA_MAX_INPUT:
                logger.warning("Total size of input (%d Mb) is more "
                               "than 100MB. Processing could take a "
                               "very long time. It is recommended to use "
                               "some other synteny backend for your data.",
                               total_size // 1024 // 1024)

            os.mkdir(work_dir)
            for block_size in self.blocks:
                block_dir = os.path.join(work_dir, str(block_size))
                perm_file = os.path.join(block_dir, "genomes_permutations.txt")
                coords_file = os.path.join(block_dir, "blocks_coords.txt")

                if not os.path.isdir(block_dir):
                    os.mkdir(block_dir)

                all_fasta = [p["fasta"] for p in recipe["genomes"].values()]
                _run_sibelia(all_fasta, block_size, block_dir)
                _postprocess_coords(chr2genome, coords_file)
                _postprocess_perms(chr2genome, perm_file)
                files[block_size] = coords_file

        return files


def _check_installation():
    return bool(utils.which(SIBELIA_EXEC))

if _check_installation():
    logger.debug("Sibelia is installed")
    SyntenyBackend.register_backend("sibelia", SibeliaBackend())
else:
    logger.debug("Sibelia is not installed")


def _get_sequence_info(recipe):
    """
    Reads fasta files and constructs a correspondence table
    between sequence name and genome name. Also returns
    the total size of input in nucleotides
    """
    chr2genome = {}
    for gen_name, gen_params in recipe["genomes"].items():
        if "fasta" not in gen_params:
            raise BackendException("FASTA file for '{0}' is not "
                                    "specified".format(gen_name))

        if not os.path.exists(gen_params["fasta"]):
            raise BackendException("Can't open '{0}'"
                                   .format(gen_params["fasta"]))

        total_size = 0
        with open(gen_params["fasta"], "r") as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    chr_name = line.strip()[1:].split(" ")[0]
                    if chr_name in chr2genome:
                        raise BackendException("Some fasta files contain "
                                               "sequences with similar names")
                    chr2genome[chr_name] = gen_name

                else:
                    total_size += len(line)

    return chr2genome, total_size


def _postprocess_perms(chr2genome, perm_file):
    """
    Converts Sibelia's permutation file to UCSC naming convention:
    genome.chromosome
    """
    new_file = perm_file + "_new"
    with open(perm_file, "r") as fin, open(new_file, "w") as fout:
        for line in fin:
            line = line.strip()
            if line.startswith(">"):
                chr_name = line[1:]
                fout.write(">{0}.{1}\n".format(chr2genome[chr_name], chr_name))
            else:
                fout.write(line + "\n")
    os.remove(perm_file)
    os.rename(new_file, perm_file)


def _postprocess_coords(chr2genome, coords_file):
    """
    Converts Sibelia's blocks_coords file to UCSC naming convention:
    genome.chromosome
    """
    new_file = coords_file + "_new"
    with open(coords_file, "r") as fin, open(new_file, "w") as fout:
        header = True
        for line in fin:
            line = line.strip()

            if header:
                if line.startswith("Seq_id"):
                    fout.write(line + "\n")
                    continue
                if line.startswith("-"):
                    fout.write(line + "\n")
                    header = False
                    continue

                chr_id, chr_size, seq_name = line.split("\t")
                fout.write("{0}\t{1}\t{2}.{3}\n".format(chr_id, chr_size,
                                                        chr2genome[seq_name],
                                                        seq_name))
            else:
                fout.write(line + "\n")

    os.remove(coords_file)
    os.rename(new_file, coords_file)


def _make_stagefile(stages, out_file):
    assert len(stages)

    with open(out_file, "w") as f:
        f.write("{0}\n".format(len(stages)))
        for stage_k, stage_d in stages:
            f.write("{0} {1}\n".format(stage_k, stage_d))


def _run_sibelia(fasta_files, block_size, out_dir):
    logger.info("Running Sibelia with block size %d", block_size)
    if not utils.which(SIBELIA_EXEC):
        raise BackendException("Sibelia is not installed")

    stagefile = os.path.join(out_dir, "stagefile.txt")
    _make_stagefile(config.vals["sibelia"], stagefile)

    devnull = open(os.devnull, "w")
    cmdline = [SIBELIA_EXEC, "--nopostprocess", "--stagefile", stagefile,
               "--minblocksize", str(block_size), "--outdir", out_dir]
    cmdline.extend(fasta_files)
    subprocess.check_call(cmdline, stdout=devnull)

    os.remove(stagefile)
    os.remove(os.path.join(out_dir, "d3_blocks_diagram.html"))
    shutil.rmtree(os.path.join(out_dir, "circos"))
