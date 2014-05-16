from __future__ import print_function
import sys

from ragout.cmaf2synteny import _make_synteny

def make_synteny(maf_file, out_dir, min_blocks_list):
    return _make_synteny(maf_file, out_dir, min_blocks_list)
