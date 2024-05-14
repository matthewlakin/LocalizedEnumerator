
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

from speciesabstract import Species_Abstract
from process import *
from strandgraph import *
import lib
import os

class FreeSpecies(Species_Abstract):

    def __init__(self, sg):
        super().__init__()
        if sg.isConnected():
            sg.__convertToCanonicalForm__()
            self.sg = sg
        else:
            errMsg = 'Tried to create a Free Species object from the following non-connected strand graph:'+os.linesep+str(self)           
            lib.error(errMsg)
            #print('ERROR: '+str(errMsg)) # Commented this out for testing purposes. Ultimately want to crash if this happens!
        # assert well-formedness check

        # assert no tethers
        for c1 in self.sg.vertex_colors:
            s = self.sg.colors_info[c1]['strand_type']
            if (s.isTethered):
                errMsg = 'Tried to create a Free Species object from the following connected strand graph that has tethers:'+os.linesep+str(self)
                #print(sg)
                lib.error(errMsg)                

    def __metric__(self):
        return self.sg.__metric__()

    def __eq__(self, other):
        assert isinstance(self, FreeSpecies)
        if isinstance(other, FreeSpecies):
            return self.__metric__() == other.__metric__()
        else:
            return False

    def __ne__(self, other):
        assert isinstance(self, FreeSpecies)
        if isinstance(other, FreeSpecies):
            return not(self.__eq__(other))
        else:
            return False

    def __lt__(self, other):
        # NB: ordering only defined between strand graphs with compatible colors!
        # NB: ordering __CURRENTLY__ only defined for connected strand graphs!
        assert isinstance(self, FreeSpecies)
        if isinstance(other, FreeSpecies):
            return self.sg.__metric__() < other.sg.__metric__()
        else: ## Should only possibly be TileSpecies?
            return False

    def __gt__(self, other):
        # NB: ordering only defined between strand graphs with compatible colors!
        # NB: ordering __CURRENTLY__ only defined for connected strand graphs!
        # assert self.compatibleColors(other)
        # assert self.isConnected()
        # assert other.isConnected()
        assert isinstance(self, FreeSpecies)
        if isinstance(other, FreeSpecies):
            return self.sg.__metric__() > other.sg.__metric__()
        else: ## Should only possibly be TileSpecies?
            return True

    # Given a connected strand graph, convert it into a Species object
    def speciesFromStrandGraph(sg):
        assert sg.isConnected()
        return FreeSpecies(sg)

    def isListOfSpecies(self, xs):
        if not isinstance(xs, list):
            return False
        for x in xs:
            if not isinstance(x, Species_Abstract):
                return False
        return True

    def connectedComponents(self):
        return self.sg.connectedComponents()

    def size(self):
        num = self.sg.numVertexes()
        return num
    
    def displayRepresentation(self):
        self.sg.displayRepresentation()    

    def new_species(self, index,  new_sg):
        return FreeSpecies(new_sg)

    def compose(self, other):
        assert isinstance(other, FreeSpecies)
        return  self.sg.compose(other.sg)

    # def speciesFromStrandGraph(self, old_strand_graph, new_strand_graph):
    #     free_species_list = []
    #     tethered_species_list = []

    #     for comp in new_strand_graph.connectedComponents():
    #         tethered_flag = False 
    #         for vertex in comp.vertex_colors:
    #             tether_info = comp.colors_info[vertex]['tether']
    #             if(tether_info[0] is not None):
    #                 tethered_flag = True               
    #         if tethered_flag:        
    #             # Add it to tethered_species list
    #             tethered_species_list.append(TileSpecies(comp))
    #         else:
    #             free_species_list.append(FreeSpecies(comp))

    #     return free_species_list + tethered_species_list 
#
# Additional helper functions for species
#

    #Check whether something is a list of species

    # #Given a , convert it into a Species object
    # def speciesfromprocess(item):
    #     sg = strandGraphComponentsFromProcess(item)
    #     assert sg.isConnected()

    #     return FreeSpecies(sg)

# Given a connected strand graph, convert it into a Species object
def speciesFromStrandGraph(sg):
    assert sg.isConnected()
    return FreeSpecies(sg)



# # Given a process, convert it into a list of free species (NB: there may be some duplicates?!)
# def speciesListFromProcess(p, domainLengthStr):
#     assert isinstance(p, Process)
#     assert p.wellFormed()
#     return [FreeSpecies(sg) for sg in connectedStrandGraphsFromProcess(p.untethered, domainLengthStr)]



# ###############################################################################################
