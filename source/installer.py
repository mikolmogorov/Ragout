#####################################
#A module for installing dependencies
#####################################

import sys, os
import subprocess
import shutil

import utils

SIBELIA_DIR="Sibelia"
SIBELIA_EXEC="Sibelia"
SIBELIA_LINK="https://github.com/bioinf/Sibelia/archive/master.tar.gz"

def install_deps(lib_dir):
    if os.path.isfile(os.path.join(lib_dir, SIBELIA_DIR, SIBELIA_EXEC)):
        sys.stdout.write("Sibelia is already installed\n")
        return True
    else:
        sys.stdout.write("Installing Sibelia\n")
        try:
            install_sibelia(lib_dir)
        except OSError:
            sys.stdout.write("Error while installing - exiting\n")
            return False

        return True


def install_sibelia(lib_dir):
    if not test_tools():
        raise OSError

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


def test_tools():
    if not utils.which("cmake"):
        sys.stdout.write("ERROR: Building Sibelia requires Cmake")
        return False

    if not utils.which("wget"):
        sys.stdout.write("ERROR: Building Sibelia requires wget")
        return False

    if not utils.which("make"):
        sys.stdout.write("ERROR: Building Sibelia requires make")
        return False

    if not utils.which("tar"):
        sys.stdout.write("ERROR: Building Sibelia requires tar")
        return False

    return True
