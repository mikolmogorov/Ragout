#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module parses Ragout configuration file
"""

from collections import namedtuple
import re
import os
import logging

from ragout.parsers.phylogeny_parser import get_leaves_names, PhyloException

logger = logging.getLogger()

class RecipeException(Exception):
    pass


def parse_ragout_recipe(filename):
    if not os.path.exists(filename):
        raise RecipeException("Can't open recipe file")

    prefix = os.path.dirname(filename)

    recipe_dict = {"genomes" : {}}
    known_params = ["tree", "target", "blocks", "maf", "fasta",
                    "circular", "draft"]
    required_params = ["tree", "target", "blocks"]

    cast_bool = ["circular", "draft"]
    cast_int_list = ["blocks"]
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
                raise RecipeException("Error parsing recipe on line {1}"
                                      .format(filename, lineno + 1))

            (obj, param_name), value = m.group(1).split("."), m.group(2)
            if param_name not in known_params:
                raise RecipeException("Unknown recipe parameter '{0}' on line {1}"
                                      .format(param_name, lineno, filename))

            #casting if necessary
            if param_name in cast_bool:
                if value in ["True", "true", "1"]:
                    value = True
                elif value in ["False", "false", "0"]:
                    value = False
                else:
                    raise RecipeException("Error parsing recipe on line "
                                          "{0}: wrong value '{1}' for bool param"
                                          .format(lineno, value))
            if param_name in cast_int_list:
                value = list(map(int, value.split(",")))
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
            raise RecipeException("Required parameter '{0}' not found in recipe"
                                  .format(param))

    genomes = None
    for param, value in recipe_dict.items():
        if param == "tree":
            try:
                genomes = get_leaves_names(value)
            except PhyloException as e:
                raise RecipeException(e)

    for g in genomes:
        recipe_dict["genomes"].setdefault(g, {})

    for g, g_params in recipe_dict["genomes"].items():
        for def_key, def_val in defaults.items():
            g_params.setdefault(def_key, def_val)

    if len(recipe_dict["blocks"]) != len(set(recipe_dict["blocks"])):
        raise RecipeException("Found similar synteny block sizes in recipe")

    if not recipe_dict["genomes"]:
        raise RecipeException("No genomes detected in recipe")

    if recipe_dict["target"] not in recipe_dict["genomes"]:
        raise RecipeException("Error parsing recipe: target genome "
                              "is not in tree")
    if "fasta" not in recipe_dict["genomes"][recipe_dict["target"]]:
        raise RecipeException("Error parsing recipe: FASTA file for "
                              "target genome is not specified")

    return recipe_dict
