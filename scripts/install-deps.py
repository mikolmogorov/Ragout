#!/usr/bin/env python

#####################################
#A script which installs dependencies
#####################################

import sys, os
import subprocess
import shutil
import argparse

import utils.utils as utils

SIBELIA_LINK = "https://github.com/bioinf/Sibelia/archive/master.tar.gz"

def install_deps(prefix):
    if utils.which("Sibelia"):
        sys.stdout.write("Sibelia is already installed\n")
        return True
    else:
        sys.stdout.write("Installing Sibelia\n")
        try:
            return install_sibelia(prefix)
        except OSError as e:
            sys.stdout.write("Error while installing - exiting\n")
            sys.stderr.write(str(e) + "\n")
            return False


def install_sibelia(prefix):
    if not test_tools():
        return False

    initial_dir = os.getcwd()

    tmp_dir = os.path.join(initial_dir, "sibelia-build")
    if os.path.isdir(tmp_dir):
        shutil.rmtree(tmp_dir)
    os.mkdir(tmp_dir)

    os.chdir(tmp_dir)
    subprocess.check_call(["wget", SIBELIA_LINK])
    subprocess.check_call(["tar", "-xf", "master.tar.gz"])

    os.chdir("Sibelia-master/build")

    srcdir = os.path.join("..", "src")
    subprocess.check_call(["cmake", srcdir, "-DONLY_SIBELIA=1",
                            "-DCMAKE_INSTALL_PREFIX=" + prefix])
    subprocess.check_call(["make"])
    subprocess.check_call(["make", "install"])


    os.chdir(initial_dir)
    shutil.rmtree(tmp_dir)

    return True


def test_tools():
    if not utils.which("cmake"):
        sys.stdout.write("ERROR: Building process requires Cmake\n")
        return False

    if not utils.which("wget"):
        sys.stdout.write("ERROR: Building process requires wget\n")
        return False

    if not utils.which("make"):
        sys.stdout.write("ERROR: Building process requires make\n")
        return False

    if not utils.which("tar"):
        sys.stdout.write("ERROR: Building process requires tar\n")
        return False

    return True


def main():
    parser = argparse.ArgumentParser(description="A helper script for"
                                                 "Sibelia installation")
    parser.add_argument("--prefix", dest="prefix",
                        help="installation prefix", default="/usr/local")
    args = parser.parse_args()
    install_deps(args.prefix)


if __name__ == "__main__":
    main()
