"""
This module provedes some functions
for debug output
"""

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

    def set_debug_dir(self, debug_dir):
        """
        Also enables debugging
        """
        self.debug_dir = debug_dir
        self.debugging = True
        if os.path.isdir(debug_dir):
            shutil.rmtree(debug_dir)
        os.mkdir(debug_dir)

    @staticmethod
    def get_instance():
        if not DebugConfig.instance:
            DebugConfig.instance = DebugConfig()
        return DebugConfig.instance
