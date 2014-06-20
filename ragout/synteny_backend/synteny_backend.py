import logging
import os
from collections import namedtuple, defaultdict

logger = logging.getLogger()


class BackendException(Exception):
    pass


class SyntenyBackend:
    backends = {}
    def __init__(self):
        pass

    def make_permutations(self, recipe, output_dir, overwrite):
        """
        Runs backend and then prepare data for futher processing
        """
        try:
            files = self.run_backend(recipe, output_dir, overwrite)
        except BackendException as e:
            logger.error(e)
            return False
        assert sorted(files.keys()) == sorted(recipe["blocks"])
        return files

    def run_backend(self, recipe, output_dir, overwrite):
        """
        Runs backend and returns a dict with permutations files
        Indexed by block sizes
        """
        return None

    @staticmethod
    def get_available_backends():
        return SyntenyBackend.backends

    @staticmethod
    def register_backend(name, instance):
        assert name not in SyntenyBackend.backends
        SyntenyBackend.backends[name] = instance
