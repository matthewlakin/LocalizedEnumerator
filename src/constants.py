
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

# Paramaters relating to DNA biophysics
DS_LENGTH = 0.34
SS_LENGTH = 0.68 # (0.68 for Chatterjee example, 0.43 for most plots from the Thubagere paper)
DSDNA_PERSISTENCE_LENGTH = 39.0
SSDNA_PERSISTENCE_LENGTH = 2.0
#HELIX_THREETURNS_LENGTH = 10.88

# Parameters related to nick angle constraints
NICKEDANGLE_UPPER_BOUND = 105
NICKED_FLAG = True

# Maximum number of sampling attempts
SAMPLING_TRIALS = 1000

# Print some/all constants
def printMainConstants(splitter=' '):
    print(f'ssDNA length per nucleotide (nm):{splitter}{DS_LENGTH}')
    print(f'dsDNA length per nucleotide (nm):{splitter}{SS_LENGTH}')
    print(f'Nick angle upper bound constraints active?{splitter}{NICKED_FLAG}')
    if NICKED_FLAG:
        print(f'Nick angle upper bound (degrees):{splitter}{NICKEDANGLE_UPPER_BOUND}')
    print(f'Max number of sampling trials:{splitter}{SAMPLING_TRIALS}')

def printAllConstants(splitter=' '):
    print(f'ssDNA persistence length (nm):{splitter}{SSDNA_PERSISTENCE_LENGTH}')
    print(f'dsDNA persistence length (nm):{splitter}{DSDNA_PERSISTENCE_LENGTH}')
    printMainConstants()
