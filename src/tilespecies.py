
from speciesabstract import Species_Abstract
from strandgraph import *
from process import *
import lib
import os


class TileSpecies(Species_Abstract):

    def __init__(self, sg_list):
        super().__init__()
        assert isinstance(sg_list, list)
        self.tiles_sg = []
        self.changedGraph = None
        for conn in sg_list:
            if conn.isConnected():
                conn.__convertToCanonicalForm__()
                self.tiles_sg.append(conn)
            else:
                errMsg = 'Tried to create a Species object from the following non-connected strand graph:'+os.linesep+str(self)
                lib.error(errMsg)
            #print('ERROR: '+str(errMsg)) # Commented this out for testing purposes. Ultimately want to crash if this happens!

        # Each strand graphs in the list contains atleast one tether and 
        for tile_species_sg in self.tiles_sg:
            tethered_flag = False
            for c in tile_species_sg.colors_info:
                if((c['tether'][0] is not None) and (c['tether'][1] is not None)):
                    tethered_flag = True
            if(not tethered_flag):
                errMsg = 'Tried to create a Tile Species object from the following connected strand graph that has no tethers:'+os.linesep+str(self)
                lib.error(errMsg)                
        self.tiles_sg.sort()

    def __metric__(self):
        assert isinstance(self, TileSpecies)
        metric_sg = []
        for tile_species in self.tiles_sg:
            metric_sg.append(tile_species.__metric__())
        return metric_sg

    def __eq__(self, other):
        assert isinstance(self, TileSpecies)
        if isinstance(other, TileSpecies):
            return (self.__metric__() == other.__metric__())
        else:
            return False
        
    def __ne__(self, other):
        assert isinstance(self, TileSpecies)
        if isinstance(other, TileSpecies):
            return not(self.__eq__(other))
        else:
            return False

    def __lt__(self, other):
        # NB: ordering only defined between strand graphs with compatible colors!
        # NB: ordering __CURRENTLY__ only defined for connected strand graphs!
        # assert self.isConnected()
        # assert other.isConnected()
        assert isinstance(self, TileSpecies)
        if isinstance(other, TileSpecies):
            return self.__metric__() < other.__metric__()
        else: ## Should only possibly be FreeSpecies?
            return True

    def __gt__(self, other):
        # NB: ordering only defined between strand graphs with compatible colors!
        # NB: ordering __CURRENTLY__ only defined for connected strand graphs!
        # assert self.compatibleColors(other)
        # assert self.isConnected()
        # assert other.isConnected()
        assert isinstance(self, TileSpecies)
        if isinstance(other, TileSpecies):
            return self.__metric__() > other.__metric__()
        else: ## Should only possibly be FreeSpecies?
            return False

    def isListOfSpecies(self, xs):
        if not isinstance(xs, list):
            return False
        for x in xs:
            if not isinstance(x, Species_Abstract):
                return False
        return True

    def removeSpeciesFromTileSpeciesList(self, sp):
        self.tiles_sg.remove(sp)

    def addSpeciesInTileSpeciesList(self, sg):
        assert sg.isConnected()
        self.tiles_sg.append(sg)

    def displayRepresentation(self):
        for sg in self.tiles_sg:
            sg.displayRepresentation()

    def connectedComponents(self):
        return self.tiles_sg

    def size(self):
        num = 0
        for sp in self.tiles_sg:
            num += sp.numVertexes()
        return num

# ###############################################################################################
