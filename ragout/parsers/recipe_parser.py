#This module parses Ragout configuration file
#############################################

from collections import namedtuple
import re
import os

RecipeParams = namedtuple("RecipeParams", ["references", "targets",
                                           "fasta", "tree", "blocks",
                                           "maf"])

class RecipeException(Exception):
    pass

#PUBLIC:
#############################################

def parse_ragout_recipe(filename):
    prefix = os.path.dirname(filename)

    ref_matcher = re.compile("REFS\s*=\s*([^\s].*)$")
    target_matcher = re.compile("TARGET\s*=\s*([^\s].*)$")
    tree_matcher = re.compile("TREE\s*=\s*([^\s].*)$")
    block_matcher = re.compile("BLOCKS\s*=\s*([^\s][\d,\s]*)$")
    fasta_matcher = re.compile("FASTA\s+(\w+)\s*=\s*([^\s]+)$")
    maf_matcher = re.compile("MAF\s*=\s*([^\s]+)$")

    tree_str = None
    block_size = None
    maf_path = None
    references = []
    targets = []
    fasta = {}

    with open(filename, "r") as f:
        for lineno, line in enumerate(f):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            m = ref_matcher.match(line)
            if m:
                references = list(map(str.strip, m.group(1).split(",")))
                continue
                #ref_id, ref_file = m.group(1), m.group(2)
                #references[ref_id] = os.path.join(prefix, ref_file)

            m = target_matcher.match(line)
            if m:
                #target_id, target_file = m.group(1), m.group(2)
                #target[target_id] = os.path.join(prefix, target_file)
                targets.append(m.group(1))
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
                    raise RecipeException("Found duplicated synteny block sizes"
                                          "while parsing {0} on line {0}"
                                          .format(filename, lineno + 1))
                continue

            m = fasta_matcher.match(line)
            if m:
                genome_name, fasta_path = m.group(1), m.group(2)
                fasta[genome_name] = os.path.join(prefix, fasta_path)
                continue

            m = maf_matcher.match(line)
            if m:
                maf_path = m.group(1)
                continue

            raise RecipeException("Error parsing {0} on line {1}"
                                  .format(filename, lineno + 1))

    if not references:
        raise RecipeException("No references specified")
    if not targets:
        raise RecipeException("No targets specified")
    if not tree_str:
        raise RecipeException("Tree is not specified")
    if not block_size:
        raise RecipeException("Synteny block sizes are not specified")

    return RecipeParams(references=references, targets=targets,
                        tree=tree_str, fasta=fasta,
                        blocks=block_size, maf=maf_path)
