
##########################################################################################
# 
# Copyright (C) 2024 Matthew Lakin, Sarika Kumar
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# 
##########################################################################################

###################################################
#
# dsdlex.py
#
# * A tokenizer for (classic) DSD species
#
# TO DO:
# ------
#
# * Increase the amount of the DSD language that can be parsed!
#    -> Will need to figure out how to expand arbitrary programs into flat lists of initial species...
#    -> Look at latest DSD manual for a grammar!
#
###################################################

import ply.lex as lex
import lib

# Lexer state for parsing comments
states = [('comment','exclusive')]

# Reserved words
reserved = { 'domain' : 'DOMAIN',
             'directive' : 'DIRECTIVE',
             'sample' : 'SAMPLE',
             'plot' : 'PLOT',
             'leak' : 'LEAK',
             'tau' : 'TAU',
             'migrate' : 'MIGRATE',
             'lengths' : 'LENGTHS',
             'def' : 'DEF',
             'new' : 'NEW',
             'true' : 'TRUE',
             'false' : 'FALSE',
             'int_of_float' : 'INT_OF_FLOAT',
             'float_of_int' : 'FLOAT_OF_INT',
             'time' : 'TIME',
             'concentration' : 'CONCENTRATION',
             'constant' : 'CONSTANT',
             'tolerance' : 'TOLERANCE',
             'sum' : 'SUM',
             'sub' : 'SUB',
             'diff' : 'DIFF',
             'div' : 'DIV',
             'scale' : 'SCALE',
             'duration' : 'DURATION',
             'points' : 'POINTS',
             'toeholds' : 'TOEHOLDS',
             'tether' : 'TETHER'
             }

# List of token names. This is required.
tokens = (['LEFTCOMMENT','RIGHTCOMMENT','LPAREN','RPAREN','LBRACE','RBRACE','LLBRACKET','RRBRACKET','LANGLE','RANGLE',
          'CARET','ASTERISK','SINGLECOLON','DOUBLECOLON','EQUALS','VERTICALBAR','UNDERSCORE','SEMICOLON',
          'ATSIGN','EXCLAMATIONMARK','INT','FLOAT', 'NUMBER','ALPHANUMERIC','SPECIALLEFTCOMMENT','SPECIALRIGHTCOMMENT', 'COMMA','PERCENT'] +
          list(reserved.values()))

# Regular expression rules for simple tokens
t_SPECIALLEFTCOMMENT  = r'\(\*\#'
t_SPECIALRIGHTCOMMENT = r'\#\*\)'
t_LPAREN              = r'\('
t_RPAREN              = r'\)'
t_LBRACE              = r'\{'
t_RBRACE              = r'\}'
t_LLBRACKET           = r'\[\['
t_RRBRACKET           = r'\]\]'
t_LANGLE              = r'\<'
t_RANGLE              = r'\>'
t_CARET               = r'\^'
t_ASTERISK            = r'\*'
t_DOUBLECOLON         = r'\:\:'
t_SINGLECOLON         = r'\:'
t_EQUALS              = r'\='
t_VERTICALBAR         = r'\|'
t_UNDERSCORE          = r'\_'
t_SEMICOLON           = r'\;'
t_ATSIGN              = r'\@'
t_EXCLAMATIONMARK     = r'\!'
t_FLOAT               = r'[0-9]+[.][0-9]*|[0-9]+[eE][+-][0-9]+|[0-9][.][0-9]*[eE][+-][0-9]+'
t_INT                 = r'[0-9]+'
t_COMMA               = r'\,'
t_PERCENT             = r'\%'

# Extract number
def t_NUMBER(t):
    r'[0-9-.]+'
    return t

# Extract reserved word, if appropriate
def t_ALPHANUMERIC(t):
    r'[A-Za-z0-9]+'
    if t.value in reserved:
        t.type = reserved.get(t.value)
    return t

# Track line numbers
def t_newline(t):
    r'\n+'
    t.lexer.lineno += len(t.value)

# Enter the 'comment' state
def t_LEFTCOMMENT(t):
    r'\(\*'
    t.lexer.level = 1        # Initial comment level
    t.lexer.begin('comment') # Enter 'comment' state

# Rules for the 'comment' state
def t_comment_LEFTCOMMENT(t):     
    r'\(\*'
    t.lexer.level +=1                
def t_comment_RIGHTCOMMENT(t):
    r'\*\)'
    t.lexer.level -=1
    if t.lexer.level == 0:
         t.lexer.begin('INITIAL')
def t_comment_newline(t):
    r'\n+'
    t.lexer.lineno += len(t.value)
def t_comment_nonspace(t):
    r'\S+'
    pass
def t_comment_error(t):
    lib.error('Illegal character \"' + str(t.value[0]) + '\" on line ' + str(t.lexer.lineno))

# A string containing ignored characters (spaces and tabs)
t_ignore  = ' \t'
t_comment_ignore = ' \t'

# Error handling rule
def t_error(t):
    lib.error('Illegal character \"' + str(t.value[0]) + '\" on line ' + str(t.lexer.lineno))

# Build the lexer
lexer = lex.lex()

# lexer.input('[[<tether (0,0) a0^*!i1%8 s*!i2%20 f^%8 s!i4%20 > | <s%20 x^%8 s*!i4%20 f^*%8 s!i2%20 a0^!i1%8>]]')

# while(True):
#     tok = lexer.token()
#     if not tok:
#         break
#     print(tok)
