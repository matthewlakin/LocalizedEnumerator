
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

# from process import *
from abc import ABC, abstractmethod
# from process import *
# from strandgraph import StrandGraph, connectedStrandGraphsFromProcess
   
###############################################################################################

#
# Species class 
#

class Species_Abstract(ABC):
    
    def __init__(self):
        super().__init__()
      
    @abstractmethod
    def __eq__(self, other):
        raise NotImplementedError

    @abstractmethod
    def __ne__(self, other):
        raise NotImplementedError

    @abstractmethod
    def __lt__(self, other):
        raise NotImplementedError

    @abstractmethod
    def __gt__(self, other):
        raise NotImplementedError

    @abstractmethod
    def __le__(self, other):
        raise NotImplementedError

    @abstractmethod
    def __ge__(self, other):
        raise NotImplementedError

    # #@abstractmethod
    # def speciesFromStrandGraph(sg):
    #     raise NotImplementedError
    
    # #@abstractmethod
    # def speciesfromprocess(self, item):
    #     raise NotImplementedError

    @abstractmethod
    def isListOfSpecies(self, specieslist):
        raise NotImplementedError
    
    @abstractmethod
    def connectedComponents(self):
        raise NotImplementedError

    @abstractmethod
    def size(self):
        raise NotImplementedError

