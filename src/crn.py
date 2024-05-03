
# 
# crn.py - a data structure for enumerated CRNs involving our species
# 

from reaction import *
#from species import isListOfSpecies # getSpeciesName, createSpeciesNames
from enumerator_geometric import *
from strandgraph import GraphvizAvailable
from freespecies import *
from tilespecies import *

###############################################################################################

class CRN(object):

    def __init__(self, species, reactions):
        self.species = species
        self.reactions = reactions
        def mkName(ctr):
            return 'sp_'+str(ctr)
        species_names = []
        for s in species:
            thisName = mkName(len(species_names))
            species_names += [(s, thisName)]
        self.species_names = species_names
        # Compress duplicate or reversible reactions within this CRN
        self.compress()
        assert self.isValid()

    # Check if this CRN is valid
    def isValid(self):
        for r in self.reactions:
            if not r.isValid():
                return False
            for s in r.listOfSpeciesInvolved():
                if s not in self.species:
                    return False
        return True

    # Go through and compress all reactions in the current CRN (combine identical ones and reversible reactions!)
    def compress(self):
        new_reactions = []
        for r in self.reactions:
            madeAChange = False
            for (sdx, s) in enumerate(new_reactions):
                x = s.tryToCombineWith(r)
                if x is not None:
                    new_reactions[sdx] = x
                    madeAChange = True
            if not madeAChange:
                new_reactions += [r]
        self.reactions = new_reactions
               
    def getSpeciesName(self, s):
        for (x,y) in self.species_names:
            if s == x:
                return y
        lib.error('In CRN.getSpeciesName, could not find species '+str(s)+' in '+str(self.species_names))
        
    def getSpecies(self, sname):
        for (x,y) in self.species_names:
            if sname == y:
                return x
        #lib.error('In CRN.getSpeciesName, could not find species '+str(sname)+' in '+str(self.species_names))
        return None


    def prettyPrintSpeciesList(self, xs):
        res = ''
        for (i,x) in enumerate(xs):
            if i > 0:
                res += ' + '
            res += self.getSpeciesName(x)
        return res

    def prettyPrintReaction(self, r): # NB: ignoring metadata here
        stringOfReactants = self.prettyPrintSpeciesList(r.reactants)
        stringOfRate = ' ->{'+str(r.fwdrate)+'} ' if r.bwdrate is None else ' {'+str(r.bwdrate)+'}<->{'+str(r.fwdrate)+'} '
        stringOfProducts = self.prettyPrintSpeciesList(r.products)
        return stringOfReactants + stringOfRate + stringOfProducts

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        def mkTitle(msg):
            return msg + os.linesep +('-'*len(msg)) + os.linesep + os.linesep
        res = mkTitle('SPECIES:')
        for s in self.species:
            res += (self.getSpeciesName(s) + ' = ' + s.printAsProcess(useNewlines=False) + os.linesep)
        res += os.linesep
        res += mkTitle('REACTIONS:')
        for r in self.reactions:
            res += (self.prettyPrintReaction(r) + os.linesep)
        return res
    
    def displayRepresentation(self):
        print('REACTIONS:')
        if self.reactions == []:
            print('None')
        else:
            for (i,r) in enumerate(self.reactions, 1):
                print('REACTION ' + str(i) + ': ' + self.prettyPrintReaction(r))
        print()
        print('KEY TO SPECIES:')
        for s in self.species:
            print(self.getSpeciesName(s)+':')
            s.displayRepresentation()


###############################################################################################
