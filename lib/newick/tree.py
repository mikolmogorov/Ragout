'''
A Python module for parsing Newick files.

Copyright (C) 2003-2008, Thomas Mailund <mailund@birc.au.dk>

This module contains the representation of trees and a parser for
creating trees from a Newick string or file. '''

import lexer
import parser


class Tree(object):
    '''
    Python representation of a tree (or rather an inner node in a tree).
    '''

    def __init__(self):
        self._edges = []
        self._leaves_cache = None
        self._identifier = None

    def add_edge(self,e):
        '''
        add_edge(e) -- add edge to sub-tree.

        Insert an edge to a new sub-tree.  The edge should be on the
        form: (st,bo,le), where st is the sub-tree, bo is the
        bootstrap value of the edge, and le is the length of the tree.
        '''
        self._edges.append(e)

        # we need to invalidate this when we add edges
        self._leaves_cache = None

    def get_edges(self):
        '''
        get_edges() -- return the list of edges to sub-trees.
        '''
        return self._edges


    def dfs_traverse(self,visitor):
        '''
        dfs_traverse(visitor) -- do a depth first traversal.

        Part of the Visitor Pattern; performs a depth first traversal,
        calling methods in visitor along the way.
        '''
        visitor.pre_visit_tree(self)
        for (n,b,l) in self._edges:
            visitor.pre_visit_edge(self,b,l,n)
            n.dfs_traverse(visitor)
            visitor.post_visit_edge(self,b,l,n)
        visitor.post_visit_tree(self)

    def get_leaves(self):
        '''
        get_leaves() --  return list of leaves in this (sub-)tree.
        '''
        if self._leaves_cache is None:
            self._leaves_cache = []
            for n,_,_ in self._edges:
                self._leaves_cache.extend(n.leaves)
        return self._leaves_cache
    def get_leaves_identifiers(self):
        """get_leaves_identifiers() --  return list of leaves' identifiers in this (sub-)tree."""
        return [l.identifier for l in self.leaves]

    @property
    def identifier(self):
        return self._identifier

    @identifier.setter
    def identifier(self, value):
        self._identifier = value

    # special functions and accessors...
    def __repr__(self):
        tree_str = '('
        sep = ''
        for (n,b,l) in self.edges:
            tree_str += sep+str(n)
            if b:
                tree_str += str(b) + ' '
            if l:
                tree_str += ' : ' + str(l)
            sep = ', '
        tree_str += ')'
        
        if self.identifier is not None:
            tree_str += self.identifier
        
        return tree_str

    edges = property(get_edges, None, None,
                     "List of edges to sub-trees.")
    leaves = property(get_leaves, None, None,
                      "List of leaves in this subtree.")
    leaves_identifiers = property(get_leaves_identifiers, None, None,
                          "List of identifiers of the leaves in this subtree.")

class Leaf(object):
    '''
    Python representation of a leaf in a tree.
    '''

    def __init__(self, identifier):
        '''
        Leaf(identifier) -- construct leaf with label identifier.
        '''
        self.identifier = identifier


    def dfs_traverse(self,visitor):
        '''
        dfs_traverse(visitor) -- do a depth first traversal.

        Part of the Visitor Pattern; calls the visit_leaf callback in visitor.
        '''
        visitor.visit_leaf(self)

    def get_leaves(self):
        '''get_leaves() --  return list of leaves in this (sub-)tree.'''
        return [self]
    def get_leaves_identifiers(self):
        """get_leaves_identifiers() --  return list of leaves' identifiers in this (sub-)tree."""
        return [self.identifier]

    def __repr__(self):
        s = self.identifier
        
        if " " in s:
            return "'"+s+"'"
        else:
            return s

    leaves = property(get_leaves, None, None,
                      "List of leaves in this subtree.")
    leaves_identifiers = property(get_leaves_identifiers, None, None,
                          "List of identifiers of the leaves in this subtree.")


class TreeVisitor(object):
    '''
    Part of the Visitor Pattern.
    '''

    def __init__(self):
        pass

    def pre_visit_tree(self,t):
        '''
        pre_visit_tree(t) -- callback called before exploring (sub-)tree t.
        '''
        pass
    def post_visit_tree(self,t):
        '''
        post_visit_tree(t) -- callback called after exploring (sub-)tree t.
        '''
        pass

    def pre_visit_edge(self,src,bootstrap,length,dst):
        '''
        pre_visit_edge(src, bo,le, dst)
        	-- callback called before exploring an edge.

        Here src is the source node and dst is the destination node,
        bo is the bootstrap support of the edge and le is the length
        of the edge.
        '''
        pass
    def post_visit_edge(self,src,bootstrap,length,dst):
        '''
        post_visit_edge(src, bo,le, dst)
        	-- callback called before exploring an edge.

        Here src is the source node and dst is the destination node,
        bo is the bootstrap support of the edge and le is the length
        of the edge.
        '''
        pass

    def visit_leaf(self,l):
        '''
        visit_leaf(l) -- callback called when exploring leaf l.
        '''
        pass


class _TreeBuilder(parser.AbstractHandler):
    def __init__(self):
        self.stack = []
        self.root = None

    def new_tree_begin(self):
        t = Tree()
        if len(self.stack) == 0:
            self.root = t
        self.stack.append(t)

    def new_tree_end(self, identifier = None):
        self.stack[-1].identifier = identifier

    def new_edge(self,bootstrap,length):
        n = self.stack.pop()
        self.stack[-1].add_edge((n,bootstrap,length))

    def new_leaf(self,l):
        if len(self.stack) == 0:        # special case of singleton tree
            self.root = l
        self.stack.append(Leaf(l))

    def get_result(self):
        return self.root


def parse_tree(input):
    '''Parse input as a Newick description of a tree and return the
    tree in a tree data structure.'''
    return parser.parse(input,_TreeBuilder())

    
def add_parent_links(tree):
    '''Extend all nodes (except for the root, of course) with a parent link.'''
    class V(TreeVisitor):
        def pre_visit_edge(self,src,b,l,dst):
            dst.parent = src
    tree.dfs_traverse(V())

def add_distance_from_root(tree):
    '''Extend all nodes with the distance (branch length) from the root'''
    tree.distance_from_root = 0.0       # 'tree' is the root...
    class V(TreeVisitor):
        def pre_visit_edge(self,src,b,l,dst):
            if l is None: l = 0
            dst.distance_from_root = src.distance_from_root + l
    tree.dfs_traverse(V())


if __name__ == '__main__':
    import unittest
    from treetest import test_suite
    unittest.TextTestRunner(verbosity=2).run(test_suite)
