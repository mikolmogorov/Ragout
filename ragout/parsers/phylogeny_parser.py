#(c) 2013-2014 by Authors
#This file is a part of Ragout program.
#Released under the BSD license (see LICENSE file)

"""
This module parses newick string and contains some helper function
to deal with trees
"""

from __future__ import absolute_import
from __future__ import division
import ragout.newick
import ragout.newick.parser
from ragout.newick.tree import Tree, Leaf

class PhyloException(Exception):
    pass


class _RagoutTreeBuilder(ragout.newick.parser.AbstractHandler):
    """
    A custom parser handler for newick library
    """
    def __init__(self):
        self.stack = []
        self.root = None

    def new_tree_begin(self):
        t = Tree()
        if len(self.stack) == 0:
            self.root = t
        self.stack.append(t)
        self.stack[-1].terminal = False

    def new_tree_end(self, identifier = None):
        self.stack[-1].identifier = identifier

    def new_edge(self,bootstrap,length):
        n = self.stack.pop()
        if length is None:
            length = 1
        self.stack[-1].add_edge((n,bootstrap,length))

    def new_leaf(self,l):
        if len(self.stack) == 0:        # special case of singleton tree
            self.root = l
        self.stack.append(Leaf(l))
        self.stack[-1].terminal = True

    def get_result(self):
        return self.root


def parse_tree(newick_str):
    tree = None
    try:
        tree = ragout.newick.parser.parse(newick_str, _RagoutTreeBuilder())
    except ragout.newick.lexer.LexerError as e:
        raise PhyloException("Error parsing tree: " + str(e))
    return tree


def get_leaves_names(newick_str):
    """
    Get names of treminal tree nodes
    """
    tree = parse_tree(newick_str)
    if tree is None:
        return None

    return [n.identifier for n in tree.get_leaves()]
