#!/usr/bin/env python

#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
A script which installs Sibelia
"""

from __future__ import print_function
import sys, os
import subprocess
import shutil
import argparse
try:
    from urllib import urlretrieve
except ImportError:
    from urllib.request import urlretrieve

SIBELIA_LINK = "https://github.com/bioinf/Sibelia/archive/master.tar.gz"
DEFAULT_PREFIX = "bin"

def install_deps(prefix):
    if which("Sibelia"):
        print("Sibelia is already installed", file=sys.stderr)
        return True
    else:
        print("Installing Sibelia", file=sys.stderr)
        try:
            return install_sibelia(prefix)
        except OSError as e:
            print("Error while installing - exiting", file=sys.stderr)
            print(e, file=sys.stderr)
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

    print("Downloading source...", file=sys.stderr)
    urlretrieve(SIBELIA_LINK, "master.tar.gz")
    subprocess.check_call(["tar", "-xf", "master.tar.gz"])

    os.chdir("Sibelia-master/build")

    srcdir = os.path.join("..", "src")
    subprocess.check_call(["cmake", srcdir, "-DONLY_SIBELIA=1",
                            "-DCMAKE_INSTALL_PREFIX=" + tmp_dir])
    subprocess.check_call(["make"])
    subprocess.check_call(["make", "install"])

    sibelia_bin_src = os.path.join(tmp_dir, "bin", "Sibelia")
    sibelia_bin_dst = os.path.join(initial_dir, prefix, "Sibelia")
    shutil.copy(sibelia_bin_src, sibelia_bin_dst)

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
    for tool in ["cmake", "make", "tar"]:
        if not which(tool):
            print("ERROR: building Sibelia requires " + tool, file=sys.stderr)
            return False
    return True


def main():
    parser = argparse.ArgumentParser(description="A helper script for "
                                                 "Sibelia installation")
    parser.add_argument("--prefix", dest="prefix",
                        help="installation prefix (default = \"{0}\")"
                        .format(DEFAULT_PREFIX),
                        default=DEFAULT_PREFIX)
    args = parser.parse_args()
    return int(not install_deps(args.prefix))


if __name__ == "__main__":
    sys.exit(main())
