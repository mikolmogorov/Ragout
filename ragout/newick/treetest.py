
import unittest
import lexer
import parser
from tree import *
from tree import _TreeBuilder

class BuilderTest(unittest.TestCase):
    ''' Test of the _TreeBuilder (and Leaf and Tree) class. '''
    
    def testTreeBuilding(self):
        ''' Test that the tree builder constructs trees correctly when
        parsed. '''
        l = lexer.Lexer("(('foo' : 0.1, 'bar' : 1.0) : 2, baz)")
        handler = _TreeBuilder()
        p = parser._Parser(l,handler)
        p.parse()
        t = handler.get_result()

        self.assertEqual(len(t.get_edges()),2)
        (t1,b1,l1), (t2,b2,l2) = t.get_edges()

        self.assertEqual(len(t1.get_edges()),2)
        self.assertEqual(l1, 2.0)
        self.assertEqual(t2.__class__, Leaf)
        self.assertEqual(l2, None)
        self.assertEqual(t.leaves_identifiers, ['foo','bar','baz'])

class TestParseTree(unittest.TestCase):
    ''' Test of the parse_tree() function. '''

    def testTreeStructure(self):
        ''' Test that a parsed tree has the right structure. '''
        t = parse_tree("(('foo' : 0.1, 'bar' : 1.0) : 2, baz)")

        self.assertEqual(len(t.get_edges()),2)
        (t1,b1,l1), (t2,b2,l2) = t.get_edges()

        self.assertEqual(len(t1.get_edges()),2)
        self.assertEqual(l1, 2.0)
        self.assertEqual(t2.__class__, Leaf)
        self.assertEqual(l2, None)
        self.assertEqual(t.leaves_identifiers, ['foo','bar','baz'])

    def testSpecialCases(self):
        ''' Test that we can parse some special cases of trees. '''

        tree = parse_tree("(B,(A,C,E),D);")
        self.assertEqual(tree.leaves_identifiers,['B','A','C','E','D'])

        tree = parse_tree("(,(,,),);")
        self.assertEqual(tree.leaves_identifiers,['']*5)

        # underscores are considered empty leaf names!
        tree = parse_tree("(_,(_,_,_),_);")
        self.assertEqual(tree.leaves_identifiers,['']*5)

        # the rest is just checking that we do not crash on this input...
        parse_tree("""
(
  ('Chimp':0.052625,
   'Human':0.042375):0.007875,
  'Gorilla':0.060125,
  ('Gibbon':0.124833,
   'Orangutan':0.0971667):0.038875
);
    """)

        parse_tree("""
(
  ('Chimp':0.052625,
   'Human':0.042375) 0.71 : 0.007875,
  'Gorilla':0.060125,
  ('Gibbon':0.124833,
   'Orangutan':0.0971667) 1.00 :0.038875
);
    """)


class TreeTest(unittest.TestCase):
    ''' Test of the Tree (and Leaf and _TreeBuilder) class. '''
    
    def testProperties(self):
        ''' Test that the tree properties lets us extract the right
        information. '''
        t = parse_tree('((A,B),C);')
        self.assertEqual(t.leaves_identifiers, ['A','B','C'])
        self.assertNotEqual(t.leaves, ['A','B','C'])

        self.assertEqual(len(t.edges), 2)
        (n1,_,_), (n2,_,_) = t.edges
        self.assertEqual(type(n1), Tree)
        self.assertEqual(type(n2), Leaf)
        self.assertEqual(n2.identifier, 'C')

class TestFunctions(unittest.TestCase):
    ''' Test of the module-level functions. '''

    def testAddParentLink(self):
        ''' Test the add_parent_links() function. '''
        t = parse_tree('((A,B),C);')
        add_parent_links(t)
        self.assertEqual([str(l.parent) for l in t.leaves],
                         ["('A', 'B')", "('A', 'B')", "(('A', 'B'), 'C')"])


    def testAddDistanceFromRoot(self):
        ''' Test the add_distance_from_root() function. '''
        t = parse_tree('((A,B),C);')
        add_distance_from_root(t)
        self.assertEqual([l.distance_from_root for l in t.leaves],[0,0,0])
        t = parse_tree('((A:2,B:3):1,C:6);')
        add_distance_from_root(t)
        self.assertEqual([l.distance_from_root for l in t.leaves],[3,4,6])


test_suite = unittest.TestSuite()
test_suite.addTest(unittest.makeSuite(BuilderTest))
test_suite.addTest(unittest.makeSuite(TestParseTree))
test_suite.addTest(unittest.makeSuite(TreeTest))
test_suite.addTest(unittest.makeSuite(TestFunctions))

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(test_suite)


#    from tree import TreeVisitor
#    def relabel(tree):
#        "Relabel the tree's leaves."
#        # visitor pattern.
#        class V(TreeVisitor):
#            def __init__(self):
#                self.count = 0
#            def visit_leaf(self,leaf):
#                leaf.identifier = str(self.count)
#                self.count += 1
#        # let visitor traverse tree
#        tree.dfs_traverse(V())
#
#    relabel(t)
#    print t
