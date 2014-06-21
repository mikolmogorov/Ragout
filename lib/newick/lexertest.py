
import unittest
import tokens
from lexer import *

class LexerTest(unittest.TestCase):
    ''' Test of the Lexer class. '''

    def test_tokens(self):
        ''' Test recognition of tokens. '''
        lexer = Lexer("()'foo' bar :0.00,;")
        lexer.read_token(tokens.LParen)
        lexer.read_token(tokens.RParen)
        identifier = lexer.read_token(tokens.ID)
        identifier = lexer.read_token(tokens.ID)
        lexer.read_token(tokens.Colon)
        n = lexer.read_token(tokens.Number)
        lexer.read_token(tokens.Comma)
        lexer.read_token(tokens.SemiColon)

    def test_identifiers(self):
        ''' Test that identifiers are correctly identified. '''
        lexer = Lexer("'foo' bar \"baz\"")
        identifier = lexer.read_token(tokens.ID)
        self.assertEqual(identifier.get_name(), 'foo')
        identifier = lexer.read_token(tokens.ID)
        self.assertEqual(identifier.get_name(), 'bar')
        identifier = lexer.read_token(tokens.ID)
        self.assertEqual(identifier.get_name(), 'baz')

    def test_numbers(self):
        ''' Test that numbers are recognized correctly. '''
        lexer = Lexer("0 0.0 0.00 1 1.0 1.11")
        n = lexer.read_token(tokens.Number)
        self.assertEqual(n.get_number(), 0.00)

test_suite = unittest.makeSuite(LexerTest)

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(test_suite)

