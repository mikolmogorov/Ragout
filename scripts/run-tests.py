#!/usr/bin/env python

#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
A script for automatic testing
"""

from __future__ import print_function
from __future__ import absolute_import
import os
import sys
import subprocess
import shutil

TESTS = {"ecoli" : {"recipe" : "examples/E.Coli/ecoli.rcp",
                    "coords" : "examples/E.Coli/mg1655.coords",
                    "max_errors" : 0,
                    "max_errors_refine" : 0,
                    "min_contigs" : 78,
                    "min_contigs_refine" : 145,
                    "max_scaffolds" : 1,
                    "outdir" : "ecoli-test",
                    "scaf_pref" : "mg1655_scaffolds"},
         #"helicobacter" : {"recipe" : "examples/H.Pylori/helicobacter.rcp",
         #                  "coords" : "examples/H.Pylori/SJM180.coords",
         #                  "max_errors" : 0,
         #                  "max_errors_refine" : 0,
         #                  "min_contigs" : 45,
         #                  "min_contigs_refine" : 126,
         #                  "max_scaffolds" : 1,
         #                  "outdir" : "helicobacter-test"},
         "cholerae" : {"recipe" : "examples/V.Cholerae/cholerae.rcp",
                       "coords" : "examples/V.Cholerae/h1.coords",
                       "max_errors" : 0,
                       "max_errors_refine" : 16,
                       "min_contigs" : 169,
                       "min_contigs_refine" : 719,
                       "max_scaffolds" : 4,
                       "outdir" : "cholerae-test",
                       "scaf_pref" : "h1_scaffolds"},
         "aureus" : {"recipe" : "examples/S.Aureus/aureus.rcp",
                     "coords" : "examples/S.Aureus/usa300.coords",
                     "max_errors" : 0,
                     "max_errors_refine" : 2,
                     "min_contigs" : 89,
                     "min_contigs_refine" : 164,
                     "max_scaffolds" : 1,
                     "outdir" : "aureus-test",
                     "scaf_pref" : "usa_scaffolds"}}

TEST_DIR = "test-dir"
RAGOUT_EXEC = "bin/ragout"
VERIFY_EXEC = os.path.join("scripts", "verify-order.py")


def test_environment():
    if not os.path.isfile(RAGOUT_EXEC):
        raise RuntimeError("File \"{0}\" was not found".format(RAGOUT_EXEC))

    if not os.path.isfile(VERIFY_EXEC):
        raise RuntimeError("File \"{0}\" was not found".format(VERIFY_EXEC))


def run_test(parameters):
    outdir = os.path.join(TEST_DIR, parameters["outdir"])
    cmd = [sys.executable, RAGOUT_EXEC, parameters["recipe"],
           "--outdir", outdir, "--debug"]
    print("Running:", " ".join(cmd), "\n")
    subprocess.check_call(cmd)

    links_simple = os.path.join(outdir, parameters["scaf_pref"] + ".links")
    links_simple_out = links_simple  + "_verify"

    #checking before refinement
    cmd = [sys.executable, VERIFY_EXEC, parameters["coords"], links_simple]
    print("Running:", " ".join(cmd), "\n")
    subprocess.check_call(cmd, stdout=open(links_simple_out, "w"))

    with open(links_simple_out, "r") as f:
        for line in f:
            if line.startswith("Total miss-ordered: "):
                value = int(line.strip()[20:])
                print("Errors:", value)
                if value > parameters["max_errors"]:
                    raise RuntimeError("Too many miss-ordered contigs")

            if line.startswith("Total contigs: "):
                value = int(line.strip()[15:])
                print("Contigs:", value)
                if value < parameters["min_contigs"]:
                    raise RuntimeError("Too few contigs")

            if line.startswith("Total scaffolds: "):
                value = int(line.strip()[17:])
                print("Scaffolds:", value)
                if value > parameters["max_scaffolds"]:
                    raise RuntimeError("Too many scaffolds")

    #checking after refinement
    cmd = [sys.executable, RAGOUT_EXEC, parameters["recipe"],
           "--outdir", outdir, "--debug", "--refine"]
    print("Running:", " ".join(cmd), "\n")
    subprocess.check_call(cmd)

    cmd = [sys.executable, VERIFY_EXEC, parameters["coords"], links_simple]
    print("Running:", " ".join(cmd), "\n")
    subprocess.check_call(cmd, stdout=open(links_simple_out, "w"))

    with open(links_simple_out, "r") as f:
        for line in f:
            if line.startswith("Total miss-ordered: "):
                value = int(line.strip()[20:])
                print("Errors:", value)
                if value > parameters["max_errors_refine"]:
                    raise RuntimeError("Too many miss-ordered contigs")

            if line.startswith("Total contigs: "):
                value = int(line.strip()[15:])
                print("Contigs:", value)
                if value < parameters["min_contigs_refine"]:
                    raise RuntimeError("Too few contigs")

            if line.startswith("Total scaffolds: "):
                value = int(line.strip()[17:])
                print("Scaffolds:", value)
                if value > parameters["max_scaffolds"]:
                    raise RuntimeError("Too many scaffolds")


def main():
    test_environment()
    if os.path.isdir(TEST_DIR):
        shutil.rmtree(TEST_DIR)
    os.mkdir(TEST_DIR)

    for name, params in TESTS.items():
        print("\n********Running test:", name, "********\n")
        run_test(params)

    print("\n********All tests were succesfully completed********")
    #shutil.rmtree(TEST_DIR)


if __name__ == "__main__":
    main()
