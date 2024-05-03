
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

    

