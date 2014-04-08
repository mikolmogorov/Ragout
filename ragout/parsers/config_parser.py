#This module parses Ragout configuration file
#############################################

from collections import namedtuple
import re
import os

ConfigParams = namedtuple("ConfigParams", ["references", "targets", "tree", "blocks"])

#PUBLIC:
#############################################

def parse_ragout_config(filename):
    prefix = os.path.dirname(filename)
    references = {}
    target = {}
    tree_str = None
    block_size = None

    ref_matcher = re.compile("REF\s+(\w+)\s*=\s*([^\s]+)$")
    target_matcher = re.compile("TARGET\s+(\w+)\s*=\s*([^\s]+)$")
    tree_matcher = re.compile("TREE\s*=\s*([^\s]+)$")
    block_matcher = re.compile("BLOCK\s*=\s*([\d,]+)$")

    for lineno, line in enumerate(open(filename, "r").read().splitlines()):
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        m = ref_matcher.match(line)
        if m:
            ref_id, ref_file = m.group(1), m.group(2)
            references[ref_id] = os.path.join(prefix, ref_file)
            continue

        m = target_matcher.match(line)
        if m:
            target_id, target_file = m.group(1), m.group(2)
            target[target_id] = os.path.join(prefix, target_file)
            continue

        m = tree_matcher.match(line)
        if m:
            tree_str = m.group(1)
            continue

        m = block_matcher.match(line)
        if m:
            sizes = m.group(1).split(",")
            block_size = list(map(int, sizes))
            if len(block_size) != len(set(block_size)):
                raise Exception("Synteny block are duplicated")
            continue

        else:
            raise Exception("Error parsing {0} on line {1}".format(filename, lineno + 1))

    return ConfigParams(references, target, tree_str, block_size)
