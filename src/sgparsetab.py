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

# sgparsetab.py
# This file is automatically generated. Do not edit.
# pylint: disable=W,C,R
_tabversion = '3.10'

_lr_method = 'LALR'

_lr_signature = 'programALPHANUMERIC ASTERISK ATSIGN CARET COMMA CONCENTRATION CONSTANT DEF DIFF DIRECTIVE DIV DOMAIN DOUBLECOLON DURATION EQUALS EXCLAMATIONMARK FALSE FLOAT FLOAT_OF_INT INT INT_OF_FLOAT LANGLE LBRACE LEAK LEFTCOMMENT LENGTHS LLBRACKET LPAREN MIGRATE NEW NUMBER PERCENT PLOT POINTS RANGLE RBRACE RIGHTCOMMENT RPAREN RRBRACKET SAMPLE SCALE SEMICOLON SINGLECOLON SPECIALLEFTCOMMENT SPECIALRIGHTCOMMENT SUB SUM TAU TETHER TIME TOEHOLDS TOLERANCE TRUE UNDERSCORE VERTICALBARbasedomain : ALPHANUMERICbasedomain : ALPHANUMERIC CARETdomain : basedomaindomain : basedomain ASTERISKdomain : basedomain EXCLAMATIONMARK ALPHANUMERICdomain : basedomain ASTERISK EXCLAMATIONMARK ALPHANUMERICnonemptydomainlist : domainnonemptydomainlist : domain nonemptydomainlistsinglestrand : LANGLE nonemptydomainlist RANGLEsinglestrand : LANGLE TETHER LPAREN NUMBER COMMA NUMBER RPAREN nonemptydomainlist RANGLEsinglestrand : LANGLE nonemptydomainlist TETHER LPAREN NUMBER COMMA NUMBER RPAREN RANGLEparallelstrand : singlestrandparallelstrand : singlestrand VERTICALBAR parallelstrandtilecontents : parallelstrandsingleobject : singlestrandsingleobject : LLBRACKET tilecontents RRBRACKETmultiobject : singleobjectmultiobject : singleobject VERTICALBAR multiobjectprocess : LPAREN multiobject RPARENprogram : process'
    
_lr_action_items = {'LPAREN':([0,15,23,],[3,24,30,]),'$end':([1,2,9,],[0,-20,-19,]),'LLBRACKET':([3,10,],[7,7,]),'LANGLE':([3,7,10,21,],[8,8,8,8,]),'RPAREN':([4,5,6,19,20,22,38,39,43,44,],[9,-17,-15,-18,-16,-9,40,41,-11,-10,]),'VERTICALBAR':([5,6,13,20,22,43,44,],[10,-15,21,-16,-9,-11,-10,]),'TETHER':([8,14,16,17,18,25,26,28,33,36,],[15,23,-7,-3,-1,-8,-4,-2,-5,-6,]),'ALPHANUMERIC':([8,16,17,18,26,27,28,32,33,36,40,],[18,18,-3,-1,-4,33,-2,36,-5,-6,18,]),'RRBRACKET':([11,12,13,22,29,43,44,],[20,-14,-12,-9,-13,-11,-10,]),'RANGLE':([14,16,17,18,25,26,28,33,36,41,42,],[22,-7,-3,-1,-8,-4,-2,-5,-6,43,44,]),'ASTERISK':([17,18,28,],[26,-1,-2,]),'EXCLAMATIONMARK':([17,18,26,28,],[27,-1,32,-2,]),'CARET':([18,],[28,]),'NUMBER':([24,30,35,37,],[31,34,38,39,]),'COMMA':([31,34,],[35,37,]),}

_lr_action = {}
for _k, _v in _lr_action_items.items():
   for _x,_y in zip(_v[0],_v[1]):
      if not _x in _lr_action:  _lr_action[_x] = {}
      _lr_action[_x][_k] = _y
del _lr_action_items

_lr_goto_items = {'program':([0,],[1,]),'process':([0,],[2,]),'multiobject':([3,10,],[4,19,]),'singleobject':([3,10,],[5,5,]),'singlestrand':([3,7,10,21,],[6,13,6,13,]),'tilecontents':([7,],[11,]),'parallelstrand':([7,21,],[12,29,]),'nonemptydomainlist':([8,16,40,],[14,25,42,]),'domain':([8,16,40,],[16,16,16,]),'basedomain':([8,16,40,],[17,17,17,]),}

_lr_goto = {}
for _k, _v in _lr_goto_items.items():
   for _x, _y in zip(_v[0], _v[1]):
       if not _x in _lr_goto: _lr_goto[_x] = {}
       _lr_goto[_x][_k] = _y
del _lr_goto_items
_lr_productions = [
  ("S' -> program","S'",1,None,None,None),
  ('basedomain -> ALPHANUMERIC','basedomain',1,'p_basedomain_long','sgparser.py',21),
  ('basedomain -> ALPHANUMERIC CARET','basedomain',2,'p_basedomain_toehold','sgparser.py',25),
  ('domain -> basedomain','domain',1,'p_domain_uncomplemented_unbound','sgparser.py',29),
  ('domain -> basedomain ASTERISK','domain',2,'p_domain_complemented_unbound','sgparser.py',33),
  ('domain -> basedomain EXCLAMATIONMARK ALPHANUMERIC','domain',3,'p_domain_uncomplemented_bound','sgparser.py',37),
  ('domain -> basedomain ASTERISK EXCLAMATIONMARK ALPHANUMERIC','domain',4,'p_domain_complemented_bound','sgparser.py',45),
  ('nonemptydomainlist -> domain','nonemptydomainlist',1,'p_nonemptydomainlist_singleton','sgparser.py',49),
  ('nonemptydomainlist -> domain nonemptydomainlist','nonemptydomainlist',2,'p_nonemptydomainlist_cons','sgparser.py',53),
  ('singlestrand -> LANGLE nonemptydomainlist RANGLE','singlestrand',3,'p_singlestrand','sgparser.py',57),
  ('singlestrand -> LANGLE TETHER LPAREN NUMBER COMMA NUMBER RPAREN nonemptydomainlist RANGLE','singlestrand',9,'p_singlestrand_tether_left','sgparser.py',61),
  ('singlestrand -> LANGLE nonemptydomainlist TETHER LPAREN NUMBER COMMA NUMBER RPAREN RANGLE','singlestrand',9,'p_singlestrand_tether_right','sgparser.py',65),
  ('parallelstrand -> singlestrand','parallelstrand',1,'p_parallelstrand_singleton','sgparser.py',69),
  ('parallelstrand -> singlestrand VERTICALBAR parallelstrand','parallelstrand',3,'p_parallelstrand_composition','sgparser.py',73),
  ('tilecontents -> parallelstrand','tilecontents',1,'p_tilecontents_composition','sgparser.py',77),
  ('singleobject -> singlestrand','singleobject',1,'p_singleobject_singelton','sgparser.py',82),
  ('singleobject -> LLBRACKET tilecontents RRBRACKET','singleobject',3,'p_singleobject_composition','sgparser.py',86),
  ('multiobject -> singleobject','multiobject',1,'p_multiobject_singelton','sgparser.py',90),
  ('multiobject -> singleobject VERTICALBAR multiobject','multiobject',3,'p_multiobject','sgparser.py',94),
  ('process -> LPAREN multiobject RPAREN','process',3,'p_process','sgparser.py',98),
  ('program -> process','program',1,'p_program','sgparser.py',102),
]
