#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module defines abstract SyntenyBackend class
"""

from __future__ import absolute_import
from __future__ import division
import logging
import os

import ragout.shared.config as config

logger = logging.getLogger()


class BackendException(Exception):
    pass


class SyntenyBackend:
    backends = {}
    def __init__(self):
        self.target_fasta = None
        self.threads = None
        self.blocks = None

    def make_permutations(self, recipe, blocks,
                          output_dir, overwrite, threads):
        """
        Runs backend and then prepare data for futher processing
        """
        self.target_fasta = recipe["genomes"][recipe["target"]].get("fasta")
        self.threads = threads
        self.blocks = blocks

        files = self.run_backend(recipe, output_dir, overwrite)
        assert sorted(files.keys()) == sorted(blocks)

        return files

    def run_backend(self, _recipe, _output_dir, _overwrite):
        """
        Runs backend and returns a dict with permutations files
        Indexed by block sizes
        """
        return None

    def get_target_fasta(self):
        """
        Returns a path to a fasta file with contigs
        """
        return self.target_fasta

    def infer_block_scale(self, recipe):
        """
        Infers synteny block scale based on target assembly size
        """
        target_fasta = recipe["genomes"][recipe["target"]].get("fasta")
        if not target_fasta or not os.path.exists(target_fasta):
            raise BackendException("Could not open target FASTA file "
                                   "or it is not specified")
        size = os.path.getsize(target_fasta)
        if size < config.vals["big_genome_threshold"]:
            return "small"
        else:
            return "large"

    @staticmethod
    def get_available_backends():
        return SyntenyBackend.backends

    @staticmethod
    def register_backend(name, instance):
        assert name not in SyntenyBackend.backends
        SyntenyBackend.backends[name] = instance
