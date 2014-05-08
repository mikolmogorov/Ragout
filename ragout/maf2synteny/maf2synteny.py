from __future__ import print_function
import sys

from ragout.cmaf2synteny import _make_synteny

def make_synteny(maf_file, out_dir, min_block):
    return _make_synteny(maf_file, out_dir, min_block)

def main():
    if len(sys.argv) != 4:
        print("Usage: maf2synteny.py maf_file out_dir min_block",
              file=sys.stderr)
        return

    maf_file = sys.argv[1]
    out_dir = sys.argv[2]
    min_block = sys.argv[3]
    make_synteny(maf_file, out_dir, int(min_block))

if __name__ == "__main__":
    main()
