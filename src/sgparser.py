
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

#
# sgparser.py - a parser for simple strand graph inputs as processes
#

######################################################################################

import ply.yacc as yacc
import lib
from domain import *
from process import *
from strand import *
from tilecontents import *

# Get the token map from the lexer. This is required.
from dsdlex import tokens

######################################################################################

def p_basedomain_long(p):
    'basedomain : ALPHANUMERIC'
    p[0] = {'domain':p[1], 'istoehold':False}

def p_basedomain_toehold(p):
    'basedomain : ALPHANUMERIC CARET'
    p[0] = {'domain':p[1], 'istoehold':True}

def p_domain_uncomplemented_unbound(p):
    'domain : basedomain'
    p[0] = Domain(p[1]['domain'], p[1]['istoehold'], False, None)

def p_domain_complemented_unbound(p):
    'domain : basedomain ASTERISK'
    p[0] = Domain(p[1]['domain'], p[1]['istoehold'], True, None)

def p_domain_uncomplemented_bound(p):
    'domain : basedomain EXCLAMATIONMARK ALPHANUMERIC'
    p[0] = Domain(p[1]['domain'], p[1]['istoehold'], False, p[3])

def p_domain_complemented_bound(p):
    'domain : basedomain ASTERISK EXCLAMATIONMARK ALPHANUMERIC'
    p[0] = Domain(p[1]['domain'], p[1]['istoehold'], True, p[4])

def p_domain_complemented_bound(p):
    'domain : basedomain ASTERISK EXCLAMATIONMARK ALPHANUMERIC'
    p[0] = Domain(p[1]['domain'], p[1]['istoehold'], True, p[4])

def p_nonemptydomainlist_singleton(p):
    'nonemptydomainlist : domain'
    p[0] = [p[1]]

def p_nonemptydomainlist_cons(p):
    'nonemptydomainlist : domain nonemptydomainlist'
    p[0] = [p[1]] + p[2]

def p_singlestrand(p):
    'singlestrand : LANGLE nonemptydomainlist RANGLE'
    p[0] = Strand(p[2], False, None, None)

def p_singlestrand_tether_left(p):
    'singlestrand : LANGLE TETHER LPAREN NUMBER COMMA NUMBER RPAREN nonemptydomainlist RANGLE'
    p[0] = Strand(p[8], True, '5prime', (p[4], p[6]))

def p_singlestrand_tether_right(p):
    'singlestrand : LANGLE nonemptydomainlist TETHER LPAREN NUMBER COMMA NUMBER RPAREN RANGLE'
    p[0] = Strand(p[2], True, '3prime', (p[5], p[7])) 

def p_parallelstrand_singleton(p):
    'parallelstrand : singlestrand'
    p[0] = [p[1]]

def p_parallelstrand_composition(p):
    'parallelstrand : singlestrand VERTICALBAR parallelstrand'
    p[0] = [p[1]] + p[3]

def p_tilecontents_composition(p):
    'tilecontents : parallelstrand'
    p[0] = TileContent(p[1])

def p_singleobject_singelton(p):
    'singleobject : singlestrand'
    p[0] = p[1] 

def p_singleobject_composition(p):
    'singleobject : LLBRACKET tilecontents RRBRACKET'
    p[0] = p[2] 

def p_multiobject_singelton(p):
    'multiobject : singleobject'
    p[0] = [p[1]]

def p_multiobject(p):
    'multiobject : singleobject VERTICALBAR multiobject'
    p[0] = [p[1]] + p[3]

def p_process(p):
    'process : LPAREN multiobject RPAREN'
    p[0] = Process(p[2])

def p_program(p):
    'program : process'
    p[0] = p[1]

# Error rule for syntax errors
def p_error(p):
    #print("Error")
    #print(p)
    if p:
         lib.error('Syntax error at token \"' + str(p.value) + '\" on line ' + str(p.lineno))
    else:
         lib.error('Syntax error at EOF')

start = 'program'

# Build the parser
parser = yacc.yacc(tabmodule="sgparsetab", errorlog=yacc.NullLogger(), debug=True)
    
# Function to enable slightly nicer use of the parser
def parse(inputtext):
    #print(inputtext)
    #print(type(parse))
    return parser.parse(inputtext)
