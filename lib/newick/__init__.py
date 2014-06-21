'''
A Python module for parsing Newick files.

Copyright (C) 2003-2008, Thomas Mailund <mailund@birc.au.dk>
'''

# convinience inclution of namespace...
from lexer  import LexerError
from parser import *
from tree   import parse_tree

if __name__ == '__main__':
    import unittest
    import lexertest
    import parsertest

    test_suite = unittest.TestSuite()
    test_suite.addTest(lexertest.test_suite)
    test_suite.addTest(parsertest.test_suite)

    unittest.TextTestRunner(verbosity=1).run(test_suite)


