
import unittest
import lexer
from parser import _Parser, AbstractHandler, parse

class ParserTest(unittest.TestCase):
    ''' Test of the _Parser class. '''
    
    def testLeafParsing(self):
        ''' Test that the parser handles leaf-nodes correctly. '''

        class LeafHandler(AbstractHandler):
            def new_leaf(self,name):
                self.name = name

        comma_lexer = lexer.Lexer(",")
        comma_handler = LeafHandler()
        comma_parser = _Parser(comma_lexer,comma_handler)
        comma_parser.parse_leaf()
        self.assertEqual(comma_handler.name,"")

        rparen_lexer = lexer.Lexer(")")
        rparen_handler = LeafHandler()
        rparen_parser = _Parser(rparen_lexer,rparen_handler)
        rparen_parser.parse_leaf()
        self.assertEqual(rparen_handler.name,"")

        number_lexer = lexer.Lexer("42")
        number_handler = LeafHandler()
        number_parser = _Parser(number_lexer,number_handler)
        number_parser.parse_leaf()
        self.assertEqual(number_handler.name,"42")

        ident_lexer = lexer.Lexer("foo")
        ident_handler = LeafHandler()
        ident_parser = _Parser(ident_lexer,ident_handler)
        ident_parser.parse_leaf()
        self.assertEqual(ident_handler.name,"foo")

        ident_lexer = lexer.Lexer("'foo'")
        ident_handler = LeafHandler()
        ident_parser = _Parser(ident_lexer,ident_handler)
        ident_parser.parse_leaf()
        self.assertEqual(ident_handler.name,"foo")

        ident_lexer = lexer.Lexer('"foo"')
        ident_handler = LeafHandler()
        ident_parser = _Parser(ident_lexer,ident_handler)
        ident_parser.parse_leaf()
        self.assertEqual(ident_handler.name,"foo")


    def testEdgeParsing(self):
        ''' Test that the parser handles edges correctly. '''

        class EdgeHandler(AbstractHandler):
            def new_edge(self,bootstrap,length):
                self.bootstrap = bootstrap
                self.length = length

        empty_edge_lexer = lexer.Lexer("dummy.leaf")
        empty_edge_handler = EdgeHandler()
        empty_edge_parser = _Parser(empty_edge_lexer,empty_edge_handler)
        empty_edge_parser.parse_edge()
        self.assertEqual(empty_edge_handler.bootstrap,None)
        self.assertEqual(empty_edge_handler.length,None)

        bootstrap_edge_lexer = lexer.Lexer("dummy.leaf 0.8")
        bootstrap_edge_handler = EdgeHandler()
        bootstrap_edge_parser = _Parser(bootstrap_edge_lexer,
                                        bootstrap_edge_handler)
        bootstrap_edge_parser.parse_edge()
        self.assertEqual(bootstrap_edge_handler.bootstrap,0.8)
        self.assertEqual(bootstrap_edge_handler.length,None)

        length_edge_lexer = lexer.Lexer("dummy.leaf : 0.5")
        length_edge_handler = EdgeHandler()
        length_edge_parser = _Parser(length_edge_lexer,
                                        length_edge_handler)
        length_edge_parser.parse_edge()
        self.assertEqual(length_edge_handler.bootstrap,None)
        self.assertEqual(length_edge_handler.length,0.5)

        full_edge_lexer = lexer.Lexer("dummy.leaf 0.8 : 0.5")
        full_edge_handler = EdgeHandler()
        full_edge_parser = _Parser(full_edge_lexer, full_edge_handler)
        full_edge_parser.parse_edge()
        self.assertEqual(full_edge_handler.bootstrap,0.8)
        self.assertEqual(full_edge_handler.length,0.5)

        # now repeating it all with a tree at the end of the edge
        empty_edge_lexer = lexer.Lexer("(dummy,tree)")
        empty_edge_handler = EdgeHandler()
        empty_edge_parser = _Parser(empty_edge_lexer,empty_edge_handler)
        empty_edge_parser.parse_edge()
        self.assertEqual(empty_edge_handler.bootstrap,None)
        self.assertEqual(empty_edge_handler.length,None)

        label_edge_lexer = lexer.Lexer("(dummy,tree)label")
        label_edge_handler = EdgeHandler()
        label_edge_parser = _Parser(label_edge_lexer,label_edge_handler)
        label_edge_parser.parse_edge()
        
        bootstrap_edge_lexer = lexer.Lexer("(dummy,tree) 0.8")
        bootstrap_edge_handler = EdgeHandler()
        bootstrap_edge_parser = _Parser(bootstrap_edge_lexer,
                                        bootstrap_edge_handler)
        bootstrap_edge_parser.parse_edge()
        self.assertEqual(bootstrap_edge_handler.bootstrap,0.8)
        self.assertEqual(bootstrap_edge_handler.length,None)

        length_edge_lexer = lexer.Lexer("(dummy,tree) : 0.5")
        length_edge_handler = EdgeHandler()
        length_edge_parser = _Parser(length_edge_lexer,
                                        length_edge_handler)
        length_edge_parser.parse_edge()
        self.assertEqual(length_edge_handler.bootstrap,None)
        self.assertEqual(length_edge_handler.length,0.5)

        full_edge_lexer = lexer.Lexer("(dummy,tree) 0.8 : 0.5")
        full_edge_handler = EdgeHandler()
        full_edge_parser = _Parser(full_edge_lexer, full_edge_handler)
        full_edge_parser.parse_edge()
        self.assertEqual(full_edge_handler.bootstrap,0.8)
        self.assertEqual(full_edge_handler.length,0.5)



    def testEdgeListParsing(self):
        ''' Test that the parser handles edge lists correctly. '''

        class EdgeListHandler(AbstractHandler):
            def __init__(self):
                super(AbstractHandler,self).__init__()

                self.bootstrap = []
                self.length = []
            def new_edge(self,bootstrap,length):
                self.bootstrap.append(bootstrap)
                self.length.append(length)

        edge_lexer = lexer.Lexer("first 1 : 0.8, second 0 : 0.5")
        edge_handler = EdgeListHandler()
        edge_parser = _Parser(edge_lexer,edge_handler)
        edge_parser.parse_edge_list()
        self.assertEqual(edge_handler.bootstrap,[1.0,0.0])
        self.assertEqual(edge_handler.length,[0.8,0.5])


    def testNodeParsing(self):
        ''' Test that the parser handles inner nodes correctly. '''

        # inner nodes are just parenthesis around an edge list, so we
        # reuse the test...
        class EdgeListHandler(AbstractHandler):
            def __init__(self):
                super(AbstractHandler,self).__init__()

                self.bootstrap = []
                self.length = []
            def new_edge(self,bootstrap,length):
                self.bootstrap.append(bootstrap)
                self.length.append(length)

        edge_lexer = lexer.Lexer("(first 1 : 0.8, second 0 : 0.5)")
        edge_handler = EdgeListHandler()
        edge_parser = _Parser(edge_lexer,edge_handler)
        edge_parser.parse_node()
        self.assertEqual(edge_handler.bootstrap,[1.0,0.0])
        self.assertEqual(edge_handler.length,[0.8,0.5])

    def testParse(self):
        ''' Test the parse() method (parsing an entire tree). '''

        class Handler(AbstractHandler):
            def __init__(self):
                super(AbstractHandler,self).__init__()
                self.edge_count = 0
                self.node_count = 0
                self.tree_end_count = 0
                self.leaf_count = 0

            def new_edge(self,bootstrap,length):
                self.edge_count += 1

            def new_tree_begin(self):
                self.node_count += 1
            def new_tree_end(self):
                self.tree_end_count += 1

            def new_leaf(self,name):
                self.leaf_count += 1


        l = lexer.Lexer("((A,B),C);")
        handler = Handler()
        parser = _Parser(l,handler)
        parser.parse_node()
        self.assertEqual(handler.node_count, handler.tree_end_count)
        self.assertEqual(handler.node_count, 2)
        self.assertEqual(handler.edge_count, 4)
        self.assertEqual(handler.leaf_count, 3)




class parseTest(unittest.TestCase):
    ''' Test of the parse() function. '''
    def testNoGetResult(self):
        ''' Test the parse() function when the handler does not have a
        get_result() method. '''

        class Handler(AbstractHandler):
            def __init__(self):
                super(AbstractHandler,self).__init__()
                self.edge_count = 0
                self.node_count = 0
                self.tree_end_count = 0
                self.leaf_count = 0

            def new_edge(self,bootstrap,length):
                self.edge_count += 1

            def new_tree_begin(self):
                self.node_count += 1
            def new_tree_end(self):
                self.tree_end_count += 1

            def new_leaf(self,name):
                self.leaf_count += 1


        handler = Handler()
        result = parse("((A,B),C);",handler)
        self.assertEqual(handler.node_count, handler.tree_end_count)
        self.assertEqual(handler.node_count, 2)
        self.assertEqual(handler.edge_count, 4)
        self.assertEqual(handler.leaf_count, 3)
        self.assertEqual(result, None)

    def testGetResult(self):
        ''' Test the parse() function when the handler does have a
        get_result() method. '''

        class Handler(AbstractHandler):
            def __init__(self):
                super(AbstractHandler,self).__init__()
                self.edge_count = 0
                self.node_count = 0
                self.tree_end_count = 0
                self.leaf_count = 0

            def new_edge(self,bootstrap,length):
                self.edge_count += 1

            def new_tree_begin(self):
                self.node_count += 1
            def new_tree_end(self):
                self.tree_end_count += 1

            def new_leaf(self,name):
                self.leaf_count += 1

            def get_result(self):
                return 42


        handler = Handler()
        result = parse("((A,B),C);",handler)
        self.assertEqual(handler.node_count, handler.tree_end_count)
        self.assertEqual(handler.node_count, 2)
        self.assertEqual(handler.edge_count, 4)
        self.assertEqual(handler.leaf_count, 3)
        self.assertEqual(result, 42)





class HandlerTest(unittest.TestCase):
    ''' Test of then handler callback protocol. '''

    class BranchLengthSum(AbstractHandler):
        def __init__(self):
            self.sum = 0.0

        def new_edge(self,bootstrap,length):
            if length:
                self.sum += length

    def test_branch_length_handler(self):
        ''' Test that the handler works, by calculating the branch length
        sum. '''
        l = lexer.Lexer("(('foo' : 0.1, 'bar' : 1.0) : 2, baz)")
        handler = HandlerTest.BranchLengthSum()
        p = _Parser(l,handler)
        p.parse()
        self.assertEqual(handler.sum, 3.1)


# FIXME: I haven't tested for correct error management!

test_suite = unittest.TestSuite()
test_suite.addTest(unittest.makeSuite(ParserTest))
test_suite.addTest(unittest.makeSuite(parseTest))
test_suite.addTest(unittest.makeSuite(HandlerTest))

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(test_suite)
    




