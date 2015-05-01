#!/usr/bin/env python2.7

#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
A script for automatic testing
"""

from __future__ import print_function
import os
import subprocess
import shutil

TESTS = {"ecoli" : {"recipe" : "examples/E.Coli/ecoli.rcp",
                    "coords" : "examples/E.Coli/mg1655.coords",
                    "max_errors" : 0,
                    "max_errors_refine" : 0,
                    "min_contigs" : 79,
                    "min_contigs_refine" : 146,
                    "max_scaffolds" : 1,
                    "outdir" : "ecoli-test"},
         "helicobacter" : {"recipe" : "examples/H.Pylori/helicobacter.rcp",
                           "coords" : "examples/H.Pylori/SJM180.coords",
                           "max_errors" : 0,
                           "max_errors_refine" : 0,
                           "min_contigs" : 45,
                           "min_contigs_refine" : 126,
                           "max_scaffolds" : 1,
                           "outdir" : "helicobacter-test"},
         "cholerae" : {"recipe" : "examples/V.Cholerae/cholerae.rcp",
                       "coords" : "examples/V.Cholerae/h1.coords",
                       "max_errors" : 0,
                       "max_errors_refine" : 13,
                       "min_contigs" : 170,
                       "min_contigs_refine" : 720,
                       "max_scaffolds" : 2,
                       "outdir" : "cholerae-test"},
         "aureus" : {"recipe" : "examples/S.Aureus/aureus.rcp",
                     "coords" : "examples/S.Aureus/usa300.coords",
                     "max_errors" : 1,
                     "max_errors_refine" : 3,
                     "min_contigs" : 100,
                     "min_contigs_refine" : 182,
                     "max_scaffolds" : 1,
                     "outdir" : "aureus-test"}}

TEST_DIR = "test-dir"
RAGOUT_EXEC = "ragout.py"
VERIFY_EXEC = os.path.join("scripts", "verify-order.py")


def test_environment():
    if not os.path.isfile(RAGOUT_EXEC):
        raise RuntimeError("File \"{0}\" was not found".format(RAGOUT_EXEC))

    if not os.path.isfile(VERIFY_EXEC):
        raise RuntimeError("File \"{0}\" was not found".format(VERIFY_EXEC))


def run_test(parameters):
    outdir = os.path.join(TEST_DIR, parameters["outdir"])
    cmd = ["python2.7", "ragout.py", parameters["recipe"],
           "--outdir", outdir, "--debug"]
    print("Running:", " ".join(cmd), "\n")
    subprocess.check_call(cmd)

    links_simple = os.path.join(outdir, "scaffolds.links")
    links_simple_out = os.path.join(outdir, "scaffolds.links_verify")
    links_refined = os.path.join(outdir, "scaffolds_refined.links")
    links_refined_out = os.path.join(outdir, "scaffolds_refined.links_verify")

    #checking before refinement
    cmd = ["python2.7", VERIFY_EXEC, parameters["coords"], links_simple]
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
    cmd = ["python2.7", VERIFY_EXEC, parameters["coords"], links_refined]
    print("Running:", " ".join(cmd), "\n")
    subprocess.check_call(cmd, stdout=open(links_refined_out, "w"))

    with open(links_refined_out, "r") as f:
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
