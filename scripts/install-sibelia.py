#!/usr/bin/env python

#####################################
#A script which installs dependencies
#####################################

from __future__ import print_function
import sys, os
import subprocess
import shutil
import argparse

SIBELIA_LINK = "https://github.com/bioinf/Sibelia/archive/master.tar.gz"

def install_deps(prefix):
    if which("Sibelia"):
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


#Mimics UNIX "which" command
def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def test_tools():
    for tool in ["cmake", "wget", "make", "tar"]:
        if not which(tool):
            print("ERROR: building Sibelia requires " + tool, file=sys.stderr)
            return False
    return True


def main():
    parser = argparse.ArgumentParser(description="A helper script for"
                                                 "Sibelia installation")
    parser.add_argument("--prefix", dest="prefix",
                        help="installation prefix", default="/usr/local")
    args = parser.parse_args()
    return int(not install_deps(args.prefix))


if __name__ == "__main__":
    sys.exit(main())
