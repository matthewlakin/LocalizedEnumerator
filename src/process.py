
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
# process.py
#

from strand import *
import lib
import os

from tilecontents import TileContent

class Process(object):

    def __init__(self, singleObjList):

        # Stores a list of tile content objects
        self.tethered = []
        # Stores a list of strands that are not tethered
        self.untethered = []
        for sobj in singleObjList:
            if(isinstance(sobj, Strand)):
                self.untethered.append(sobj)
            elif(isinstance(sobj, TileContent)):
                self.tethered.append(sobj)
            else:
                assert False

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return self.compactString()

    # (Maybe) compact string representation
    def compactString(self, useNewlines=False):
        output = '( '
        for (idx,s) in enumerate(self.untethered):
            if idx > 0:
                output += ' | '
            output += str(s)
            if idx == len(self.untethered) - 1:
                if useNewlines:
                    output += os.linesep
        for i, tc in enumerate(self.tethered) :
            for (idx,s) in enumerate(tc.strands):
                if (idx > 0 or len(self.untethered) > 0):
                    output += ' | '
                if(idx == 0):
                    output += '[[ '
                output += str(s) # This works now, but produces ugly output as this is actually a strand graph!
            if i == len(self.tethered) - 1:
                output += ' ]]'
                if useNewlines:
                    output += os.linesep                  
        output += ' )'
        return output

    # We assume that all processes P are well-formed in that each
    # bond i in P appears exactly twice and is shared between complementary domains.
    def wellFormed(self):
        bonds_dict = {} 
        for s in self.untethered:
            for d in s.domains:
                if d.bond is not None:
                    bondname = d.bond
                    if bondname not in bonds_dict:
                        bonds_dict[bondname] = [d]
                    else:
                        bonds_dict[bondname] += [d]

        for tile_content in self.tethered:
            for s in tile_content.tethered_strands:
                for d in s.domains:
                    if d.bond is not None:
                        bondname = d.bond
                        if bondname not in bonds_dict:
                            bonds_dict[bondname] = [d]
                        else:
                            bonds_dict[bondname] += [d]

        for (bondname, domains) in bonds_dict.items():
            if (len(domains) == 2) and (domains[0].wellFormedBondTo(domains[1])):
                pass
            else:
                print("bond names are not correct!")
                return False

        
        ## Return false if strands has tether and is not inside double square brackett [[ ]] or not tethered on a tile.
        for s in self.untethered: 
            if(s.isTethered or s.tether_orientation is not None  or s.tether_coord is not None):
                print("Strand is in untethered list, but tether information is provided.")
                return False
            else:
                pass

        unprocessed_list = [] #self.tethered
        for tile_content in self.tethered:
            for s in tile_content.tethered_strands:
                unprocessed_list.append(s)
        attempted_list = {}
        tethered_bond = {}

        while(len(unprocessed_list) > 0):
            s = unprocessed_list.pop(0)
            if(s.isTethered and s.tether_orientation is not None and s.tether_coord is not None):                
                for d in s.domains:
                    if d.bond is not None:
                        bondname = d.bond
                        if bondname not in tethered_bond:
                            tethered_bond[bondname] = True
            else:
                flag_tethered = False
                for d in s.domains:
                    if d.bond is not None:
                        bondname = d.bond
                        if bondname in tethered_bond:
                            flag_tethered = True
                            break
                if(flag_tethered):
                    for d in s.domains:
                        if d.bond is not None:
                            bondname = d.bond
                            tethered_bond[bondname] = True       
                else:
                    unprocessed_list.append(s)
                    if(str(s) in attempted_list.keys()):
                        if(attempted_list[str(s)] == 10):
                            print("strand is not tethered: "+str(s))
                            return False
                        else:
                            attempted_list[str(s)] += 1
                    else:
                        attempted_list[str(s)] = 1
        # print("Well formed True!!!")
        return True                    

        # for tile_content in self.tethered:
        #     for s in tile_content.tethered_strands:
        #         if (s.isTethered):
        #             pass
        #         elif(not self.isStrandTethered(s, [], tile_content.tethered_strands)):
        #                 return False
        #         else:
        #             pass
        
        # for tile_content in self.tethered:
        #     tether_coord = []
        #     for s in tile_content.tethered_strands:
        #         if (s.isTethered):
        #             if (s.tether_coord in tether_coord):
        #                 return False
        #             else:
        #                 tether_coord.append(s.tether_coord)
        # return True

    # def isStrandTethered(self, s, strand_list, strands):
    #     if (s in strand_list):
    #         return False
    #     for d in s.domains:
    #         if (d.bond is not None):
    #             d_comp = d.getComplement()
    #             st = self.getStrand(d_comp, strands)
    #             if (st is not None):
    #                 if (st.isTethered):
    #                     return True
    #                 else:
    #                     if(self.isStrandTethered(st, strand_list, strands)):
    #                         return True
    #     strand_list.append(s)
    #     #return False

    # def getStrand(self, d, strands):
    #     for s in strands:  
    #         for dom in s.domains:
    #             if(d == dom):
    #                 return s
    #     return None
         
