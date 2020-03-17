#!/usr/bin/env python

#(c) 2019 by Authors
#This file is a part of the Flye package.
#Released under the BSD license (see LICENSE file)

"""
Runs simple toy test
"""


from __future__ import print_function

import os
import sys
import subprocess
import shutil
from distutils.spawn import find_executable


def test_toy():
    if not find_executable("ragout"):
        sys.exit("ragout is not installed!")

    print("Running toy test:\n")
    script_dir = os.path.dirname(os.path.realpath(__file__))
    recipe_file = os.path.join(script_dir, "data", "ecoli.rcp")
    out_dir = "raout_toy_test"
    subprocess.check_call(["ragout", recipe_file, "-o", out_dir, "--refine"])
    shutil.rmtree(out_dir)
    print("\nTEST SUCCESSFUL")


def main():
    test_toy()
    return 0


if __name__ == "__main__":
    sys.exit(main())
