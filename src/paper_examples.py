
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

import sgparser
from speciesabstract import *
from enumerator_geometric import *
from strandgraph import *
import math
from constraintchecker_sampling import *
import sys
from timeit import default_timer as timer
from constants import *

defaultSeed = 11

def mkEnumeratorGeometric(seed):
    return ReactionEnumerator_Geometric({'name':'adjacent_detailed',
                                         'debug': False,
                                         'enumerationMode':'detailed',
                                         'maxComplexSize': math.inf,
                                         'threeWayMode':'adjacent',
                                         'unbindingMode':'adjacent',
                                         'rate' : {'bind': 0.003, 'unbind' : 0.1, 'migrate' : 1.0, 'displace' : 1.0},
                                         'constraintChecker': ConstraintChecker_Sampling(seed=seed)})

def process_input(s, domainLengthStr, enumeratorGeometric, verbose=True):
    if verbose:
        print('PARAMETERS:')
        printMainConstants()
    p = sgparser.parse(s)
    species_list = speciesFromProcess(p, domainLengthStr)
    if verbose:
        print()
        print("INITIAL SPECIES:")
        for index, species in enumerate(species_list):
            print(f'sp_{index}:')
            species.displayRepresentation()
        print()
    crn = enumeratorGeometric.enumerateReactions(species_list) # Actually do the enumeration!
    if verbose:
        print('Found '+str(len(crn.species))+' species and '+str(len(crn.reactions))+' reactions in total.')
        print()
    return crn

########################################################################

# Entry point for Chatterjee (2017) localized strand displacement examples
def chatterjee_circuit(dist_between_hairpins=10.88, seed=defaultSeed, verbose=True):
    domainLengthStr = 'toeholdDomain spcr1 length 5 toeholdDomain spcr2 length 5 longDomain s length 12 toeholdDomain a0 length 6 toeholdDomain f length 6 toeholdDomain x length 6 longDomain y length 12'
    s = '( <s a0^> | [[<tether(0,0) spcr1 a0^* s*!i1 f^ s!i1> | <tether('+ str(dist_between_hairpins) +',0) spcr2 x^* s*!i3 y^ s!i3> ]] | <s!i2 x^ s*!i2 f^*> | <s*!i4 y^*> | <s!i4> )'
    start_time = timer()
    if verbose:
        print(f'INPUT STRINGS:\n{domainLengthStr}\n{s}\n')
    enumeratorGeometric = mkEnumeratorGeometric(seed)
    crn = process_input(s, domainLengthStr, enumeratorGeometric, verbose=verbose)
    if verbose:
        crn.displayRepresentation()
    end_time = timer()
    elapsed_time = end_time - start_time
    if verbose:
        print('Time taken to enumerate reactions for settings '+enumeratorGeometric.settings['name']+': '+str(elapsed_time)+' seconds')
        print()
    return crn

# Entry point for Thubagere (2017) cargo-sorting robot examples
def enumerate_robot_cargo(s, domainLengthStr, robot, cargo, seed=defaultSeed, verbose=True):
    start_time = timer()
    print(f'INPUT STRINGS:\n{domainLengthStr}\n{s}\n')
    enumeratorGeometric = mkEnumeratorGeometric(seed)
    crn = process_input(s, domainLengthStr, seed, verbose=verbose)
    robot_sg = strandGraphFromProcess(sgparser.parse(robot), domainLengthStr)
    cargo_sg = strandGraphFromProcess(sgparser.parse(cargo), domainLengthStr)
    robot_strand = robot_sg.colors_info[0]['strand_type']
    cargo_strand = cargo_sg.colors_info[0]['strand_type']
    crn.displayModifiedRepresentation(robot_strand, cargo_strand)
    end_time = timer()
    elapsed_time = end_time - start_time
    print('Time taken to enumerate reactions for settings '+enumeratorGeometric.settings['name']+': '+str(elapsed_time)+' seconds')
    print()
    return crn

