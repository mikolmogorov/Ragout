#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module parses newick string and contains some helper function
to deal with trees
"""

import newick

class PhyloException(Exception):
    pass


def parse_tree(newick_str):
    tree = None
    try:
        tree = newick.parse_tree(newick_str)
    except newick.lexer.LexerError as e:
        raise PhyloException("Error parsing tree")
    return tree


def get_leaves_names(newick_str):
    """
    Get names of treminal tree nodes
    """
    tree = parse_tree(newick_str)
    if tree is None:
        return None

    return list(map(lambda n: n.identifier, tree.get_leaves()))


def is_terminal(tree_node):
    """
    A kind of hack to identify leaf
    """
    return tree_node.get_leaves() == [tree_node]
