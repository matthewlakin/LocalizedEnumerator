
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
# enumerator_abstract.py - abstract superclass for defining reaction enumerators
#

from abc import ABC, abstractmethod
from process import Process
from speciesabstract import *
from freespecies import *
from tilespecies import *
import sgparser

class ReactionEnumerator_Abstract(ABC):

    def __init__(self):
        super().__init__()
        
    #
    # ABSTRACT METHOD:
    # Given a list of distinct species, enumerate all possible reactions that could occur.
    #
    @abstractmethod
    def enumerateReactions(self, species_list):
        pass

    def enumerateReactionsForProcess(self, p):
        assert type(p) == Process
        return self.enumerateReactions(speciesListFromProcess(p))

    def enumerateReactionsForInputString(self, s):
        return self.enumerateReactionsForProcess(sgparser.parse(s))
