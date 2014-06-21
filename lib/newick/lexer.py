'''
A Python module for parsing Newick files.

Copyright (C) 2003-2008, Thomas Mailund <mailund@birc.au.dk>

This module contains the functionality for lexical analysis.  You
rarely need to access it in your own modules and should probably only
use the parser or tree modules, and from those most likely only the
objects loaded into the mail package.  '''

import tokens
import re

_patterns = [
    (tokens.Number, 	re.compile(r'\s*(-?\d+(\.\d+)?([eE][+-]?\d+)?)\s*(?=[,:(); \t\n])')),
    (tokens.ID, 	re.compile(r"\s*((\"[^\"]+\")|('[^']+')|(\w[^,:(); \t\n]*|_)+)\s*")),
    (tokens.Colon, 	re.compile(r'\s*(:)\s*')),
    (tokens.SemiColon, 	re.compile(r'\s*(;)\s*')),
    (tokens.Comma, 	re.compile(r'\s*(,)\s*')),
    (tokens.LParen, 	re.compile(r'\s*(\()\s*')),
    (tokens.RParen, 	re.compile(r'\s*(\))\s*'))
    ]

class LexerError(Exception):
    '''Exception thrown if the lexer encounters an error.'''
    def __init__(self,err):
        self.err = err

    def __str__(self):
        return "LexerError: "+self.err


class Lexer(object):
    '''Lexicographical analysis of a Newick tree.'''

    def __init__(self, input):
        self.input = input
        self.next_token = None

    def remaining(self):
        ''' The remaining input stream, i.e. the stream that hasn't been split
        into tokens. '''
        result = None
        if self.next_token:
            result = str(self.next_token)+" "+self.input
        else:
            result = self.input
        result.strip()
        return result

    def peek_next_token(self):
        ''' Returns the next token in the input, without deleting it
        from the input stream. '''
        if self.next_token:
            return self.next_token
        else:
            for (cons, p) in _patterns:
                m = re.match(p,self.input)
                if m:
                    self.next_token = cons(self.input[m.start():m.end()])
                    self.input = self.input[m.end():]
                    return self.next_token
            # no match, either end of string or lex-error
            if self.input:
                raise LexerError("Unknown token at "+self.input[:10]+"...")
            else:
                return None

    def get_next_token(self):
        ''' Returns (and delete) the next token from the input
        stream. '''
        token = self.peek_next_token()
        self.next_token = None
        return token

    def read_token(self,token_class):
        ''' Read a token of the specified class, or raise an exception
        if the next token is not of the given class. '''
        token = self.get_next_token()
        if token.__class__ != token_class:
            raise LexerError("expected "+str(token_class)+
                             " but received "+str(token.__class__)+
                             " at "+self.input[:10]+"...")
        else:
            return token

    def peek_token(self,token_class):
        ''' Checks whether the next token is of the specified class. '''
        token = self.peek_next_token()
        return token.__class__ == token_class



if __name__ == '__main__':
    import unittest
    from lexertest import test_suite
    unittest.TextTestRunner(verbosity=2).run(test_suite)
