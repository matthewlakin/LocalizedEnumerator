
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
# crn_modified.py - a data structure for enumerated CRNs involving our species
# 

from reaction import *
from crn import CRN
from enumerator_geometric import *
from strandgraph import GraphvizAvailable
from tilespecies import *

###############################################################################################

class CRN_Modified(CRN):

    def __init__(self, species, reactions):
        super().__init__(species, reactions)

    def modifiedPrettyPrintReaction(self, r, robot, cargo): 
        robot_cargo_info_reac = self.getRobotandCargoInfo(r.reactants, robot, cargo)
        robot_cargo_info_prod = self.getRobotandCargoInfo(r.products, robot, cargo)
        stringOfReactants = self.prettyPrintSpeciesList(r.reactants)
        stringOfRate = ' ->{'+str(r.fwdrate)+'} ' if r.bwdrate is None else ' {'+str(r.bwdrate)+'}<->{'+str(r.fwdrate)+'} '
        stringOfProducts = self.prettyPrintSpeciesList(r.products)
        res = stringOfReactants
        if (len(robot_cargo_info_reac) > 0):
            if(robot is not None):
                res +=  ' @robot:' + str(robot_cargo_info_reac[0]) 
            if(cargo is not None):
                res += ' @cargo:' + str(robot_cargo_info_reac[1])
        res +=  stringOfRate + stringOfProducts
        if (len(robot_cargo_info_prod) > 0):
            if(robot is not None):
                res +=  ' @robot:' + str(robot_cargo_info_prod[0]) 
            if(cargo is not None):
                res += ' @cargo:' + str(robot_cargo_info_prod[1])        
        return res

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
            #res += (self.modifiedPrettyPrintReaction(r) + os.linesep)
            res += (self.prettyPrintReaction(r) + os.linesep)
        return res
    
    def displayModifiedRepresentation(self, robot, cargo):
        self.displayGraphicalRepresentation(robot, cargo)
        print('REACTIONS:')
        if self.reactions == []:
            print('None')
        else:
            for (i,r) in enumerate(self.reactions, 1):
                print('REACTION ' + str(i) + ': ' + self.modifiedPrettyPrintReaction(r, robot, cargo))
        print()
        print('KEY TO SPECIES:')
        for s in self.species:
            print(self.getSpeciesName(s)+':')
            s.displayRepresentation()

    def getRobotandCargoInfo(self, species_list, robot, cargo):
        #In Thubagre's paper, robot strand is always bounded to tracks. 
        #Here, we are only interested in the tether information of tracks which are tilespecies.
        robot_tether_info = []
        cargo_tether_infor = []
        for sp in species_list:
            if (isinstance(sp, TileSpecies)):
                robot_info, cargo_info = self.getRobotandCargoInfo_fromSpecies(sp, robot, cargo)
                robot_tether_info += robot_info
                cargo_tether_infor += cargo_info
        return robot_tether_info, cargo_tether_infor

    def getRobotandCargoInfo_fromSpecies(self, sp, robot, cargo):
        assert (isinstance(sp, TileSpecies))
        robot_tether_coord = []
        cargo_tether_coord = []
        for tile_sg in sp.tiles_sg:
            # tile_sg is a connected strand graph
            assert tile_sg.isConnected()
            vs = tile_sg.vertex_colors
            flag_cargo_strand = False
            flag_robot_strand = False
            tether_coord = []
            for v in vs:
                tether = tile_sg.colors_info[v]['tether']
                strand_type = tile_sg.colors_info[v]['strand_type'] 
                if (robot is not None and strand_type == robot):
                    flag_robot_strand = True
                elif (cargo is not None and strand_type == cargo):
                    flag_cargo_strand = True
                else:
                    pass
                if (tether[0] is not None and tether[1] is not None):
                        if(tether[1] not in tether_coord):
                            tether_coord += [tether[1]]
            if(flag_cargo_strand):
                cargo_tether_coord += tether_coord
            if(flag_robot_strand):
                robot_tether_coord += tether_coord
        return (robot_tether_coord, cargo_tether_coord)

    def makeGraphicalRepresentation(self, robot, cargo):
        if not GraphvizAvailable:
            lib.error('Graphviz module not available!')
        
        d = graphviz.Digraph(engine='dot',strict=True)
        #d.attr(colorscheme='X11')
        #d.attr(rankdir='LR')
        #d.attr(overlap_scaling='15.0')

        for sp in self.species:
            str_sp = self.getSpeciesName(sp)
            (robot_loc, cargo_loc) = self.getRobotandCargoInfo_fromSpecies(sp, robot, cargo)
            if(len(robot_loc) != 0 and len(cargo_loc) != 0):
                label = str_sp + "\\n robot: "+str(robot_loc) + '\\n cargo: '+str(cargo_loc)
                d.node(str_sp, label, color='black', fontcolor='black', fontname='Arial', fontsize="10" ,shape='rectangle')
            elif(len(cargo_loc) != 0):
                label = str_sp + " \\n cargo: "+str(cargo_loc)
                d.node(str_sp, label, color='black', fontcolor='black', fontname='Arial', fontsize="10" ,shape='rectangle')
            elif(len(robot_loc) != 0):
                label = str_sp + " \\n robot: "+str(robot_loc)
                d.node(str_sp,  label, color='black', fontcolor='black', fontname='Arial', fontsize="10" ,shape='rectangle')
            else:
                d.node(str_sp, str_sp, color='black', fontcolor='black', fontname='Arial', fontsize="10" ,shape='rectangle')

        flag_reversible = False
        for r1 in self.reactions:
            if r1.bwdrate is None :
                flag_reversible = False
            else :
                flag_reversible = True
            for reac in r1.reactants:
                str_reac = self.getSpeciesName(reac)
                for product in r1.products:
                    str_prod = self.getSpeciesName(product)
                    
                    if(flag_reversible):
                        d.edge(str_reac, str_prod, color='black') 
                        d.edge(str_prod, str_reac, color='black')
                    else:
                        d.edge(str_reac, str_prod, color='red') 
        return d

    # Display graphical representation of this strand graph.
    def displayGraphicalRepresentation(self, robot, cargo):
        try:
            testJupyter = display
            useGraphviz = True
        except NameError: # Should jump back to here if not in Jupyter, as display won't be defined by default. NB: any other NameError will be picked up too, however...
            useGraphviz = False
        if useGraphviz:
            display(self.makeGraphicalRepresentation(robot, cargo))
        else:
            print(self)

###############################################################################################
