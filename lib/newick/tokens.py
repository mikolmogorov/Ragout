'''
A Python module for parsing Newick files.

Copyright (C) 2003-2008, Thomas Mailund <mailund@birc.au.dk>

This module contains the tokens used in the parser. '''

class Token(object):
    def __init__(self, str):
        self.str = str

    def __repr__(self):
        return 'T"'+self.str+'"'

class LParen(Token):
    pass

class RParen(Token):
    pass

class ID(Token):
    def __init__(self, identifier):
        identifier = identifier.strip()
        if identifier[0] in ("'",'"'):
            identifier = identifier[1:-1]
        self.identifier = identifier

    def get_name(self):
        return self.identifier

    def __repr__(self):
        return 'ID"'+self.identifier+'"'

class Colon(Token):
    pass

class SemiColon(Token):
    pass

class Comma(Token):
    pass

class Number(Token):
    def __init__(self, number):
        if "." in number:
            self.number = float(number)
        else:
            self.number = int(number)

    def get_number(self):
        return self.number

    def __repr__(self):
        return 'NUMBER"'+str(self.number)+'"'
