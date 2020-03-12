"""
This module provedes some functions
for debug output
"""

from __future__ import absolute_import
from __future__ import division
import os
import shutil


class DebugConfig():
    """
    Singleton providing global debug configuration
    """
    instance = None

    def __init__(self):
        self.debug_dir = None
        self.debugging = False

    def set_debugging(self, debugging):
        self.debugging = debugging

    def set_debug_dir(self, debug_dir):
        if not self.debugging:
            return
        self.debug_dir = debug_dir
        if not os.path.isdir(debug_dir):
            os.mkdir(debug_dir)

    def clear_debug_dir(self):
        if not self.debugging:
            return
        if os.path.isdir(self.debug_dir):
            shutil.rmtree(self.debug_dir)
            os.mkdir(self.debug_dir)

    @staticmethod
    def get_instance():
        if not DebugConfig.instance:
            DebugConfig.instance = DebugConfig()
        return DebugConfig.instance
