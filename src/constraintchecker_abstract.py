from abc import ABC, abstractmethod

class ConstraintChecker_Abstract(ABC):

    def __init__(self):
        super().__init__()
        
    #
    # ABSTRACT METHOD:
    
    def isPlausible(self, sg):
        pass
