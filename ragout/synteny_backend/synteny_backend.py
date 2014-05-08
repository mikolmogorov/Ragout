import logging
import os
from collections import namedtuple, defaultdict

logger = logging.getLogger()

#PUBLIC:
####################################################

class BackendException(Exception):
    pass


class SyntenyBackend:
    backends = {}
    def __init__(self):
        pass

    #runs backend and then prepare data for futher processing
    def make_permutations(self, config, output_dir, overwrite):
        try:
            files = self.run_backend(config, output_dir, overwrite)
        except BackendException as e:
            logger.debug(e)
            return False
        assert sorted(files.keys()) == sorted(config.blocks)
        return files

    #runs backend and returns a dict with permutations files
    #indexed by block sizes
    def run_backend(self, config, output_dir, overwrite):
        return None

    @staticmethod
    def get_available_backends():
        return SyntenyBackend.backends

    @staticmethod
    def register_backend(name, instance):
        assert name not in SyntenyBackend.backends
        SyntenyBackend.backends[name] = instance
