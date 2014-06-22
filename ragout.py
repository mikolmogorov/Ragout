#!/usr/bin/env python2.7

#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This script does all the necessary preparations
and invokes Ragout
"""

import os
import sys

LIB_DIR = "lib"

#Check Python version
if sys.version_info[:2] != (2, 7):
    print("Error: Ragout requires Python version 2.7 ({0}.{1} detected)."
          .format(sys.version_info[0], sys.version_info[1]))
    sys.exit(-1)

#Setting executable paths
ragout_root = os.path.dirname(os.path.realpath(__file__))
lib_absolute = os.path.join(ragout_root, LIB_DIR)
sys.path.extend([ragout_root, lib_absolute])
os.environ["PATH"] += os.pathsep + lib_absolute

#Ragout entry point
from ragout.main import main
sys.exit(main())
