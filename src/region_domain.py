
########################################################################

from constants import *

########################################################################

# Basic representation of individual domains: just ssDNA vs dsDNA and length in nucleotides
class RegionDomain:
    def __init__(self, isDS, lengthNT):
        self.isDS = isDS
        self.lengthNT = lengthNT
    
    def __repr__(self):
        return ('DS' if self.isDS else 'SS') + str(self.lengthNT)

    def __str__(self):
        return self.__repr__()

    def maxLength(self):
        return self.lengthNT * (DS_LENGTH if self.isDS else SS_LENGTH)

    def persistenceLength(self):
        return (DSDNA_PERSISTENCE_LENGTH
                if self.isDS
                else SSDNA_PERSISTENCE_LENGTH)

########################################################################
