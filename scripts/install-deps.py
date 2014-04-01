#!/usr/bin/env python2

#####################################
#A script which installs dependencies
#####################################

import sys, os
import subprocess
import shutil

parent_dir = os.path.dirname(os.path.abspath(os.path.dirname(__file__)))
import utils.utils as utils

LIB_DIR = os.path.join(parent_dir, "lib")
SIBELIA_DIR = "Sibelia"
SIBELIA_EXEC = "Sibelia"
SIBELIA_LINK = "https://github.com/bioinf/Sibelia/archive/master.tar.gz"

def install_deps(lib_dir):
    if os.path.isfile(os.path.join(lib_dir, SIBELIA_DIR, SIBELIA_EXEC)):
        sys.stdout.write("Sibelia is already installed\n")
        return True
    else:
        sys.stdout.write("Installing Sibelia\n")
        try:
            return install_sibelia(lib_dir)
        except OSError as e:
            sys.stdout.write("Error while installing - exiting\n")
            sys.stderr.write(str(e) + "\n")
            return False


def install_sibelia(lib_dir):
    if not test_tools():
        return False

    initial_dir = os.getcwd()

    sibelia_dir = os.path.abspath(os.path.join(lib_dir, SIBELIA_DIR))
    if not os.path.isdir(lib_dir):
        os.mkdir(lib_dir)
    if not os.path.isdir(sibelia_dir):
        os.mkdir(sibelia_dir)

    tmp_dir = os.path.join(sibelia_dir, "build-tmp")
    if os.path.isdir(tmp_dir):
        shutil.rmtree(tmp_dir)
    os.mkdir(tmp_dir)

    os.chdir(tmp_dir)
    subprocess.check_call(["wget", SIBELIA_LINK])
    subprocess.check_call(["tar", "-xf", "master.tar.gz"])

    install_dir = os.path.abspath("Sibelia-install")
    os.mkdir(install_dir)
    os.chdir("Sibelia-master/build")

    srcdir = os.path.join("..", "src")
    subprocess.check_call(["cmake", srcdir, "-DONLY_SIBELIA=1",
                            "-DCMAKE_INSTALL_PREFIX=" + install_dir])
    subprocess.check_call(["make"])
    subprocess.check_call(["make", "install"])

    shutil.copy(os.path.join(install_dir, "bin", SIBELIA_EXEC),
                os.path.join(sibelia_dir, SIBELIA_EXEC))

    os.chdir(initial_dir)
    shutil.rmtree(tmp_dir)

    return True


def test_tools():
    if not utils.which("cmake"):
        sys.stdout.write("ERROR: Building Sibelia requires Cmake\n")
        return False

    if not utils.which("wget"):
        sys.stdout.write("ERROR: Building Sibelia requires wget\n")
        return False

    if not utils.which("make"):
        sys.stdout.write("ERROR: Building Sibelia requires make\n")
        return False

    if not utils.which("tar"):
        sys.stdout.write("ERROR: Building Sibelia requires tar\n")
        return False

    return True


if __name__ == "__main__":
    install_deps(LIB_DIR)
