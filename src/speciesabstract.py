
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
# species.py - specific subset of strand graphs that are connected, and therefore represent a single species.
#            - currently implemented as a subclass of StrandGraph!
#            - NB: we could do away with this if can generalize graph enumeration to non-connected strand graphs...
#

from process import *
from abc import ABC, abstractmethod
from process import *
# from strandgraph import StrandGraph, connectedStrandGraphsFromProcess
   
###############################################################################################

#
# Species class 
#

class Species_Abstract(ABC):
    
    def __init__(self):
        super().__init__()
      
    # ABSTRACT METHOD:
    # Given a list of distinct species,

    # Given a connected strand graph, convert it into a Species object
    def speciesFromStrandGraph(sg):
        pass
    
    def speciesfromprocess(self, item):
        pass

    def isListOfSpecies(self, specieslist):
        pass

    def connectedComponents(self):
        pass

    def size(self):
        pass

    

