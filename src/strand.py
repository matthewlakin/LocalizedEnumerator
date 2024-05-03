
from domain import *

class Strand(object):

    def __init__(self, domains, isTethered, tether_orientation, tether_coord):
        assert isinstance(domains, list)
        for d in domains:
            assert isinstance(d, Domain)
        self.domains = domains
        self.isTethered = isTethered
        self.tether_orientation = tether_orientation
        self.tether_coord = (float(tether_coord[0]), float(tether_coord[1])) if tether_coord is not None else None

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        output = '<'
        if (str(self.tether_orientation) == '5prime'):
            output += ' '  +str(self.tether_orientation) + ' ' +str(self.tether_coord) + ' '

        for (idx, d) in enumerate(self.domains):
            if idx != 0:
                output += ' '
            output += str(d)
        
        if (str(self.tether_orientation) == '3prime'):
            output += ' '  +str(self.tether_orientation) + ' ' +str(self.tether_coord)
        output += '>'
        return output

    def __metric__(self):
        return (self.domains, self.tether_orientation, self.tether_coord)

    def __eq__(self, other):
        assert isinstance(other, Strand)
        return (self.__metric__() == other.__metric__())

    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        assert isinstance(other, Strand)
        return (self.__metric__() < other.__metric__())

    def __gt__(self, other):
        assert isinstance(other, Strand)
        return (self.__metric__() > other.__metric__())

    def strandType(self):
        return Strand([d.stripBond() for d in self.domains], self.isTethered, self.tether_orientation, self.tether_coord)

    def copyStrand(self):
        return Strand([d.copyDomain() for d in self.domains], self.isTethered, self.tether_orientation, self.tether_coord)

    def modifyDomain(self, idx, newDomain):
        assert 0 <= idx <= len(self.domains) - 1
        new_domains = [d.copyDomain() for d in self.domains]
        new_domains[idx] = newDomain
        return Strand(new_domains, self.isTethered, self.tether_orientation, self.tether_coord)
