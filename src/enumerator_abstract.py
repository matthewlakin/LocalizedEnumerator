
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
