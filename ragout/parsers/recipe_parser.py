#This module parses Ragout configuration file
#############################################

from collections import namedtuple
import re
import os
import logging
from Bio import Phylo
try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO

logger = logging.getLogger()
#RecipeParams = namedtuple("RecipeParams", ["genomes", "targets",
#                                           "tree", "blocks",
#                                           "maf"])
#Genome

class RecipeException(Exception):
    pass

#PUBLIC:
#############################################

def parse_ragout_recipe(filename):
    if not os.path.exists(filename):
        raise RecipeException("Can't open recipe file")

    prefix = os.path.dirname(filename)

    recipe_dict = {"genomes" : {}}
    required_params = ["tree", "target", "blocks"]
    known_params = ["tree", "target", "blocks", "maf", "fasta",
                    "circular", "draft"]

    cast_bool = ["circular", "draft"]
    cast_int_list = ["blocks"]
    cast_str_list = ["target"]
    fix_path = ["fasta", "maf"]

    defaults = {"circular" : False,
                "draft" : False}

    param_matcher = re.compile("([^\s]+)\s*=\s*([^\s].*)$")
    with open(filename, "r") as f:
        for lineno, line in enumerate(f):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            m = param_matcher.match(line)
            if not m or not "." in m.group(1):
                raise RecipeException("Error parsing {0} on line {1}"
                                      .format(filename, lineno + 1))

            (obj, param_name), value = m.group(1).split("."), m.group(2)
            if param_name not in known_params:
                raise RecipeException("Unknown parameter '{0}' on line {1} "
                                      "while parsing {2}"
                                      .format(param_name, lineno, filename))

            #casting if necessary
            if param_name in cast_bool:
                value = bool(value)
            if param_name in cast_int_list:
                value = list(map(int, value.split(",")))
            if param_name in cast_str_list:
                value = value.split(",")
            if param_name in fix_path:
                value = os.path.join(prefix, value)

            if obj == "":
                recipe_dict[param_name] = value
            elif obj == "*":
                defaults[param_name] = value
            else:
                recipe_dict["genomes"].setdefault(obj, {})[param_name] = value

    for param in required_params:
        if param not in recipe_dict:
            raise RecipeException("Required parameter '{0}' not found in config"
                                  .format(param))

    genomes = None
    for param, value in recipe_dict.items():
        if param == "tree":
            tree = Phylo.read(StringIO(value), "newick")
            genomes = list(map(lambda n: n.name, tree.get_terminals()))

    for g in genomes:
        recipe_dict["genomes"].setdefault(g, {})

    for g, g_params in recipe_dict["genomes"].items():
        for def_key, def_val in defaults.items():
            g_params.setdefault(def_key, def_val)

    if len(recipe_dict["blocks"]) != len(set(recipe_dict["blocks"])):
        raise RecipeException("Found duplicated synteny block sizes"
                              "while parsing {0}".format(filename))

    if not recipe_dict["genomes"]:
        raise RecipeException("No genomes in terminal nodes of the tree")

    for t in recipe_dict["target"]:
        if t not in recipe_dict["genomes"]:
            raise RecipeException("Target genome is not in tree")
        if "fasta" not in recipe_dict["genomes"][t]:
            logger.warning("FASTA file for target genome is not specified")

    return recipe_dict
