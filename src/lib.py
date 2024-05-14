
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

######################################################
#
# lib.py - general utility functions
#
######################################################

import random
import string
import inspect
import sys
import subprocess
import os.path
import datetime
import re
import itertools
import math

# Detect whether running Python 2 or Python 3
inPython3 = sys.version_info >= (3,0)

# Get current date codes
def getDateCode():
    return datetime.datetime.now().strftime("%A %d %B %Y").rstrip()
    ## return datetime.datetime.now().strftime("%m-%d-%Y")
def getTimeCode():
    return datetime.datetime.now().strftime("%I:%M %p %Z").rstrip()
    ## return datetime.datetime.now().strftime("%H:%M:%S")

# Print a float to specified number of significant figures (but don't round it!)
def fstr(f, dp=1):
    assert type(dp) is int
    assert dp >= 1
    formatStr = '{0:.'+str(dp)+'f}'
    return formatStr.format(f)

# Round to significant figures
def round_sig(x, sig=2):
    if x == 0.0:
        return 0.0
    else:
        return round(x, -int(math.floor(math.log10(math.fabs(x)))) + (sig - 1)) 

# Are the list elements all distinct?
def distinct(xs):
    temp = []
    for x in xs:
        if x in temp:
            return False
        else:
            temp.append(x)
    return True

# Do two lists have any elements in common?
def overlap(xs,ys):
    for x in xs:
        for y in ys:
            if x == y:
                return True
    return False

# Identity function
def identity(x):
    return x

# Produce a string from a list, using the given function and given separator
def string_of_list(xs, f, sep):
    res = ''
    for (i,x) in enumerate(xs):
        if i == 0:
            res += f(x)
        else:
            res += (sep + f(x))
    return res

# Flatten a list of lists
def flatten(xss):
    return [x for xs in xss for x in xs]

# All permutations of a list, as a list of lists
def allPermsAsLists(xs):
    return [list(p) for p in itertools.permutations(xs)]

# Generate a random string of letters, with length 10 by default
def random_letters(length=10):
    prefix = ''
    for i in range(0,length):
        allLetters = string.ascii_letters if inPython3 else string.letters
        prefix += random.choice(allLetters)
    return prefix
def random_prefix(length=10):
    return random_letters(length=length)

# Generate a randomly-named folder in the current directory, keep checking until a new name is found
def makeRandomFolder(prefix=None, separator='_', length=10, withinFolder='.'):
    while True:
        suffix = random_letters(length=length)
        if prefix is not None:
            name = prefix + separator + suffix
        else:
            name = suffix
        if not os.path.exists(withinFolder + os.sep + name):
            os.mkdir(withinFolder + os.sep + name)
            return name

# Concatenate a list of strings
def concat(xs):
    return ''.join(xs)

# Reverse a string
def reverse(str):
    return str[::-1]

# Display an error message and quit
def error(x):
    print('ERROR: ' + x)
    sys.exit(1)

# Signal an internal error and quit
def internalError():
    file = inspect.getfile(inspect.currentframe().f_back)
    lineno = inspect.currentframe().f_back.f_lineno
    print('INTERNAL ERROR in file ' + str(file) + ', line ' + str(lineno))
    sys.exit(1)

# Try to execute a system call
def trySystemCall(cmd):
    if subprocess.call(cmd, shell=True) != 0:
        error('Failed to execute command "' + cmd + '"')
def trySystemCall_Quiet(cmd):
    if subprocess.call(cmd, shell=True, stdout=subprocess.DEVNULL) != 0:
        error('Failed to execute command "' + cmd + '"')    
def trySystemCall_Silent(cmd):
    if subprocess.call(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) != 0:
        error('Failed to execute command "' + cmd + '"')

# Try to execute a system call, and capture the output
def trySystemCallWithOutput(cmd):
    try:
        res = subprocess.check_output(cmd, shell=True)
        return res.decode('utf-8').split(os.linesep)
    except CalledProcessError:
        error('Failed to execute command "' + cmd + '"')

# Debug printing
debugMode = False
def activateDebugMode():
    global debugMode
    debugMode = True
def disableDebugMode():
    global debugMod
    debugMode = False
def debug(x):
    if debugMode:
        print(x)

# Higher-order list for_all function
def for_all(xs, p):
    for x in xs:
        if not p(x):
            return False
    return True

# Higher-order list exists function
def exists(xs, p):
    for x in xs:
        if p(x):
            return True
    return False

# Does a string contain a valid, well-nested dotparen string?
# Check that all characters are valid, parens are well-nested and all strands are connected
def valid_dotparen(s):
    assert isinstance(s, str)
    if s == '':
        return False
    nestLevel = 0
    if exists(s, lambda c: c not in ['.', '(', ')', '+']):
        return False
    numStrandBreaks = 0
    for c in s:
        if c == '.':
            pass
        elif c == '+':
            numStrandBreaks += 1
            if nestLevel == 0:
                return False
        elif c == '(':
            nestLevel += 1
        elif c == ')':
            if nestLevel > 0:
                nestLevel -= 1
            else:
                return False
        else:
            return False
    if nestLevel != 0:
        return False
    return True

# Number of nucleotides in dotparen structure
def nucleotidesInDotparen(s):
    assert valid_dotparen(s)
    return len(s.replace('+',''))

# Does a string contain a valid, well-nested DU+ string?
# Assumes all tokens are whitespace-separated (except "Dn(" tokens)...
# Checks that all characters are valid, parens are well-nested and all strands are connected
duplusUnboundRE = re.compile(r'^U[1-9][0-9]*$')
duplusOpenDuplexRE = re.compile(r'^D[1-9][0-9]*\($')
def valid_duplus(s):
    assert isinstance(s, str)
    if s == '':
        return False
    nestLevel = 0
    numStrandBreaks = 0
    tokens = s.split()
    for token in tokens:
        if duplusUnboundRE.match(token) is not None:
            pass
        elif token == '+':
            numStrandBreaks += 1
            if nestLevel == 0:
                return False
        elif duplusOpenDuplexRE.match(token) is not None:
            nestLevel += 1
        elif token == ')':
            if nestLevel > 0:
                nestLevel -= 1
            else:
                return False
        else:
            return False
    if nestLevel != 0:
        return False
    return True

# Number of nucleotides in DU+ structure
def nucleotidesInDUPlus(s):
    assert valid_duplus(s)
    tokens = s.split()
    numNucleotides = 0
    nestings = []
    for token in tokens:
        if duplusUnboundRE.match(token) is not None:
            numNucleotides += int(token[1:])
        elif duplusOpenDuplexRE.match(token) is not None:
            num = int(token[1:-1])
            nestings += [num]
            numNucleotides += num
        elif token == ')':
            num = nestings.pop(-1)
            numNucleotides += num
        elif token == '+':
            pass
        else:
            assert False
    assert nestings == []
    return numNucleotides

def removeStrandBreaks(dp):
    return dp.replace('+','')

def numBasesInDotParen(dp):
    assert valid_dotparen(dp)
    return len(removeStrandBreaks(dp))

def numBasesInDecoratedDotParen(decdp):
    assert valid_decorated(decdp)
    return len(decdp)

# Find all indexes of a character in a string
def findChar(s, ch):
    return [i for i, ltr in enumerate(s) if ltr == ch]

# Decorate dotparen to annotate with the base pairs etc
def decorate_dotparen(s):
    assert valid_dotparen(s)
    ## print('Initial dotparen: ' + s)
    strandBreaksAfter = [i-jdx-1 for (jdx,i) in enumerate(findChar(s,'+'))]
    s = removeStrandBreaks(s)
    ## print('Dotparen with strand breaks removed: ' + s)
    ## print('strandBreaksAfter: ' + str(strandBreaksAfter))
    openingParens = []
    basePairsDict = {}
    for (idx,c) in enumerate(s):
        if c == '(':
            openingParens.append(idx)
        elif c == ')':
            basePairsDict[openingParens[-1]] = idx
            basePairsDict[idx] = openingParens[-1]
            openingParens.pop()
    assert openingParens == []
    ## print('basePairsDict: ' + str(basePairsDict))
    results = []
    for (idx,c) in enumerate(s):
        isFinalBase = idx == len(s)-1
        nextBackboneIs3PrimeEnd = isFinalBase or (idx in strandBreaksAfter)
        if c == '.':
            boundToIndex = None
        elif c == '(' or c == ')':
            boundToIndex = basePairsDict[idx]
        else:
            internalError()
        this_result = {'index':idx, 'boundToIndex':boundToIndex, 'nextBackboneIs3PrimeEnd':nextBackboneIs3PrimeEnd}
        results.append(this_result)
    return results

# Is argument a valid decorated dot-paren structure?
def valid_decorated(decdp):
    for d in decdp:
        if sorted(list(d.keys())) != sorted(['index','boundToIndex','nextBackboneIs3PrimeEnd']):
            return False
        if d['index'] < 0 or d['nextBackboneIs3PrimeEnd'] not in [True,False]:
            return False
        if d['boundToIndex'] == None:
            continue
        elif d['boundToIndex'] < 0:
            return False
    return True

# Does a string contain only valid nucleotides (and is non-empty)?
def valid_nucleotides(s):
    return for_all(s, lambda c: c in ['A','C','G','T']) and len(s) > 0

# Does a string contain only valid nucleotide codes (and is non-empty)?
def valid_codes(s):
    return for_all(s, lambda c: c in ['A','C','G','T','N','R','Y','M','K','S','W','V','H','B','D','U']) and len(s) > 0

# Do two characters represent complementary bases, using just Watson-Crick pairs?
wc_dict = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
def complementary_bases(c1,c2):
    return (c2 == wc_dict[c1])

# Do two characters represent complementary bases, if we also allow G-T wobble pairs?
wc_wobble_dict = {'A':['T'], 'C':['G'], 'G':['C','T'], 'T':['A','G']}
def complementary_bases_wobble(c1,c2):
    return (c2 in wc_wobble_dict[c1])

# Translate codes into list of possible nucleotides
codes_dict = {'A':['A'], 'C':['C'], 'G':['G'], 'T':['T'], 'N':['A','C','G','T'], 'R':['A','G'], 'Y':['C','T'], 'M':['A','C'],
              'K':['G','T'], 'S':['C','G'], 'W':['A','T'], 'V':['A','C','G'], 'H':['A','C','T'], 'B':['C','G','T'], 'D':['A','G','T'],
              'U':['T']}
def translate_code(c):
    if c in codes_dict:
        return codes_dict[c]
    else:
        error('Nucleic acid code ' + c + ' not recognised')

# Convert (DNA) nucleotide sequence into corresponding RNA sequence string suitable for IDT ordering
def get_IDT_RNA_string(s):
    assert valid_nucleotides(s)
    idt_rna_translation_dict = {'A':'rA', 'C':'rC', 'G':'rG', 'T':'rU'}
    idt_individual_rna_bases = [idt_rna_translation_dict[b] for b in s]
    return ''.join(idt_individual_rna_bases)

# Get a single random base, optionally excluding a subset of the bases
def getRandomBase(exclude=[]):
    allowed_nucleotides = ['A','C','G','T']
    assert valid_nucleotides(allowed_nucleotides)
    assert valid_nucleotides(exclude)
    for x in exclude:
        if x in allowed_nucleotides:
            allowed_nucleotides.remove(x)
    if allowed_nucleotides == []:
        assert False
    else:
        return random.choice(allowed_nucleotides)

# Return a permutation of the supplied nucleotide sequence
def scrambleNucleotides(seq):
    assert valid_nucleotides(seq)
    return random.sample(seq, k=len(seq))

# Mutate a randomly chosen base in a strand
def mutateRandomBase(strand):
    assert valid_nucleotides(strand)
    idx = random.randint(0, len(strand)-1)
    currentBase = strand[idx]
    newBase = getRandomBase(exclude=[currentBase])
    newStrand = strand[0:idx]+newBase+strand[idx+1:]
    return newStrand

# Compute G-C percentage in single strand
def computeGCPercentage(strand):
    num_gc = (strand.count('G') + strand.count('C'))
    num_bases = len(strand)
    return (float(num_gc) / float(num_bases)) * 100.0

# Compute total G-C percentage in list of strands
def computeTotalGCPercentage(strands):
    num_gc = 0
    num_bases = 0
    for strand in strands:
        assert valid_nucleotides(strand)
        num_gc += (strand.count('G') + strand.count('C'))
        num_bases += len(strand)
    return (float(num_gc) / float(num_bases)) * 100.0

# Check whether a sequence matches the supplied template (Must be same length!)!
def sequenceMatches(seq,template):
    assert valid_nucleotides(seq)
    assert valid_codes(template)
    assert len(seq) == len(template)
    for (base,code) in list(zip(seq,template)):
        if base not in translate_code(code):
            return False
    return True

# Return all subsequence windows of given length
def windows(x, wsize, wstep=1):
    assert wsize <= len(x)
    assert wstep > 0
    return list([x[i:i+wsize] for i in range(0,len(x)-wsize+1,wstep)])

# Return all subsequence windows of given length, along with their starting indexes
def indexed_windows(x, wsize, wstep=1):
    assert wsize <= len(x)
    assert wstep > 0
    return list([(i,x[i:i+wsize]) for i in range(0,len(x)-wsize+1,wstep)])

# Check whether a sequence contains a subsequence matching the supplied template!
def anySubsequenceMatches(seq,template):
    assert valid_nucleotides(seq)
    assert valid_codes(template)
    if len(seq) >= len(template):
        for subseq in windows(seq, len(template)):
            if sequenceMatches(subseq,template):
                return True
        return False
    else:
        return False

# Return all subsequences from a sequence that match the supplied template!
def getAllSubsequenceMatches(seq,template):
    assert valid_nucleotides(seq)
    assert valid_codes(template)
    matches = []
    if len(seq) >= len(template):
        for subseq in windows(seq, len(template)):
            if sequenceMatches(subseq,template):
                matches.append(subseq)
    return matches

# Check whether a sequence contains no subsequence matching the supplied template!
def noSubsequenceMatches(seq,template):
    return not anySubsequenceMatches(seq,template)

# Check whether a sequence contains no subsequences matching any of the supplied templates!
def noSubsequenceMatchesAnyTemplate(seq, templates):
    for template in templates:
        if anySubsequenceMatches(seq,template):
            return False
    return True

# Generate the reversed (Watson-Crick) complement of a sequence
def rev_comp(s):
    assert valid_nucleotides(s)
    s = reverse(s) # Reverse s
    out = ''
    for c in s:
        if c in wc_dict:
            out += wc_dict[c]
        else:
            raise Exception('invalid string')
    return out

# Change one base in a sequence, as specified
def change_one_base(s, idx, newbase):
    assert valid_nucleotides(s)
    assert len(s) > 0
    assert valid_nucleotides(newbase)
    assert len(newbase) == 1
    assert 0 <= idx < len(s)
    if idx == 0:
        res = newbase + s[1:]
    elif idx == len(s)-1:
        res = s[:-1] + newbase
    else:
        res = s[:idx] + newbase + s[idx+1:]
    assert len(res) == len(s)
    return res

# Return a random sequence conforming to the supplied template
def random_sequence(template):
    assert valid_codes(template)
    if valid_nucleotides(template):
        # Special case when there is no freedom to vary sequence
        return template
    if len(template) == 1:
        return random.choice(translate_code(template[0]))
    else:
        return random.choice(translate_code(template[0])) + random_sequence(template[1:])

# Return a random sequence conforming to the supplied template, avoiding certain subsequences
def random_sequence_avoiding(template, avoids):
    while True:
        seq = random_sequence(template)
        if noSubsequenceMatchesAnyTemplate(seq, avoids):
            return seq

# Return all sequences matching supplied template
def all_sequences(template):
    if len(template) == 1:
        return translate_code(template[0])
    else:
        shorter_seqs = all_sequences(template[1:])
        these_seqs = []
        for s in shorter_seqs:
            for b in translate_code(template[0]):
                these_seqs.append(b + s)
        return these_seqs

# Return string representation of a sequence
def sequence_string(s, rev=False, margin=0, ends=True):
    if len(s) < 1:
        error('Empty sequence passed to sequence_string')
    x = ' '*margin
    rng = list(range(0,len(s)))
    if rev:
        if ends:
            x += '3\'-'
        rng = reverse(rng)
    else:
        if ends:
            x += '5\'-'
    for i in rng:
        x += s[i]
    if ends:
        if rev:
            x += '-5\''
        else:
            x += '-3\''
    return x

# Generate random concrete domains according to the supplied domain templates.
def random_concrete_domains(dtemplates):
    cdomains = {}
    for (d,t) in list(dtemplates.items()):
        cdomains[d] = random_sequence(t)
    return cdomains

# Read a file, and return each line in a list
def read_file(fname):
    lines = []
    with open(fname, 'r') as f:
        for line in f:
            lines.append(line)
    return lines

# Write lines to a file
def write_file(fname, lines, quiet=False):
    with open(fname, 'w') as f:
        for line in lines:
            f.write(line + os.linesep)
    if not quiet:
        print('Wrote output to '+fname)

# Read and tokenize a file, and return a list of lists of strings.
def tokenize(fname):
    lines = []
    for line in read_file(fname):
        lines.append(line.split())
    return lines

# Trim comments out from a line, expressed as a list of strings.
# We consider a comment to start from any token whose initial character is '#', and extend to the rest of the line
def stripHashComments(line):
    res = []
    for token in line:
        if token == '':
            continue
        elif token[0] == '#':
            break
        else:
            res.append(token)
    return res

# See whether token is a float
def isFloat(x):
    try:
        a = float(x)
        return True
    except ValueError:
        return False

# See whether token is an int (or a float that represents an actual integer)
def isInt(x):
    try:
        a = float(x)
        b = int(a)
    except ValueError:
        return False
    else:
        return a == b

# See whether a token is an int OR a float
def isNumber(token):
    return isFloat(token) or isInt(token)

# Get a float from a token, if possible
def getFloat(token):
    assert isFloat(token)
    return float(token)

# Does a token represent None?
def isNone(token):
    return token == 'None'

# Is a token either float or None?
def isFloatOrNone(token):
    return isFloat(token) or isNone(token)

# Get a float from a token, or None
def getFloatOrNone(token):
    assert isFloatOrNone(token)
    if isNone(token):
        return None
    else:
        return float(token)

# Get an int from a token that represents an int
def getInt(token):
    assert isInt(token)
    return int(token)

# Is a token either int or None?
def isIntOrNone(token):
    return isInt(token) or isNone(token)

# Get an int from a token, or None
def getIntOrNone(token):
    assert isIntOrNone(token)
    if isNone(token):
        return None
    else:
        return int(token)

# Check whether a token is a boolean
def isBool(token):
    return token in ['true','false','True','False']

# Get a boolean from a token, if possible
def getBool(token):
    assert isBool(token)
    if token in ['true','True']:
        return True
    elif token in ['false','False']:
        return False
    else:
        assert False

# Line is empty
def lineIsEmpty(line):
    return line == []

# Try to match a specific token
def matchToken(token, expected):
    return token in expected

# See whether token is equals
def isEquals(token):
    return matchToken(token,['='])

# Strip specified string (stripme, 2nd arg) from beginning of xs (1st arg), if it occurs there.
# This is different from xs.lstrip(arg) because '1233456'.lstrip('123') returns '456',
# whereas lstripString('1233456', '123') returns '3456'.
def lstripString(xs, stripme):
    if stripme == '':
        return xs
    else:
        if xs.startswith(stripme):
            return xs[len(stripme):]
        else:
            return xs

# Strip specified string (stripme, 2nd arg) from end of xs (1st arg), if it occurs there.
# This is different from xs.rstrip(arg) because '1233456'.rstrip('3456') returns '12',
# whereas rstripString('1233456', '3456') returns '123'.
def rstripString(xs, stripme):
    if stripme == '':
        return xs
    else:
        if xs.endswith(stripme):
            return xs[:-len(stripme)]
        else:
            return xs

# Discard any final comma
def discardFinalComma(token):
    return rstripString(token, ',')

# Split the list xs into n sub-lists of roughly equal size
def divide(xs, n):
    assert n <= len(xs)
    csize = len(xs) / n
    acc = []
    for i in range(0,n):
        acc.append((xs[(csize*i):(csize*(i+1))]))
    if csize*n < len(xs):
        for x in xs[(csize*n):]:
            acc[-1].append(x)
    return acc

# Basic Hamming distance function
# def hammingDistBasic(xs,ys):
#     assert len(xs) == len(ys)
#     h = 0
#     for (x,y) in list(zip(xs,ys)):
#         if x != y:
#             h += 1
#     return h

# Compute Hamming distance, return None if there is a matching substring of length greater than threshold
def hammingDistAdvanced(xs,ys,maxMatchLength=None):
    assert len(xs) == len(ys)
    assert maxMatchLength == None or maxMatchLength > 0
    longestMatchLength = 0
    currentMatchLength = 0
    h = 0
    for (x,y) in list(zip(xs,ys)):
        if x != y:
            h += 1
            if currentMatchLength > longestMatchLength:
                longestMatchLength = currentMatchLength
                currentMatchLength = 0
        else:
            currentMatchLength += 1
    if currentMatchLength > longestMatchLength:
        longestMatchLength = currentMatchLength
    ## print('Hamming distance was ' + str(h))
    ## print('Longest match length was ' + str(longestMatchLength))
    if maxMatchLength == None:
        return h
    else:
        if longestMatchLength > maxMatchLength:
            ## print('PROBLEM: longest match length was greater than allowed threshold (' + str(maxMatchLength) + ')')
            return None
        else:
            ## print('Longest match length was less than or equal to allowed threshold (' + str(maxMatchLength) + ')')
            return h

# Calculates the Levenshtein distance of 2 strings (http://code.activestate.com/recipes/576874/)
def levenshtein(s1, s2):
    l1 = len(s1)
    l2 = len(s2)
    matrix = [list(range(l1 + 1))] * (l2 + 1)
    for zz in range(l2 + 1):
        matrix[zz] = list(range(zz,zz + l1 + 1))
    for zz in range(0,l2):
        for sz in range(0,l1):
            if s1[sz] == s2[zz]:
                matrix[zz+1][sz+1] = min(matrix[zz+1][sz] + 1, matrix[zz][sz+1] + 1, matrix[zz][sz])
            else:
                matrix[zz+1][sz+1] = min(matrix[zz+1][sz] + 1, matrix[zz][sz+1] + 1, matrix[zz][sz] + 1)
    return matrix[l2][l1]

# Return a list of pairs of neighboring elements from a list (of even length)
def getNeighboringPairs(l):
    assert len(l) % 2 == 0
    return [(l[i],l[i+1]) for i in range(0,len(l),2)]

# Get input from the command-line (different for Python 2 vs Python 3)
def get_input(prompt):
    assert isinstance(prompt, str)
    if inPython3:
        return input(prompt)
    else:
        return raw_input(prompt)

# Get a yes response from the user
# Returns True if user enters something starting with 'y' or 'Y' at the prompt, and False otherwise
def get_yes_input(prompt):
    response = get_input(prompt)
    return len(response) > 0 and response[0] in ['y','Y']

# Find (first) key with given value in dict
def findKeyWithGivenValue(d,v):
    for (k,vprime) in list(d.items()):
        if vprime == v:
            return k
    return None

# Split a list on all elements that satisfy supplied predicate (with option to keep / discard the separators)
def splitListOnPredicate(xs, p, keepSeparators=False, keepEmpties=True):
    allSubparts = []
    thisSubpart = []
    for x in xs:
        if p(x):
            if (keepEmpties) or ((not keepEmpties) and thisSubpart != []):
                allSubparts.append(list(thisSubpart))
            if keepSeparators:
                allSubparts.append([x])
            thisSubpart = []
        else:
            thisSubpart.append(x)
    if (keepEmpties) or ((not keepEmpties) and thisSubpart != []):
        allSubparts.append(list(thisSubpart))
    return allSubparts

# Split into pieces of given lengths (must be correct length)
def splitIntoLengths(xs, lens):
    assert len(xs) == sum(lens)
    iacc = 0
    acc = []
    for l in lens:
        acc += [xs[iacc:iacc+l]]
        iacc += l
    return acc

# Return list of pairs of pairwise neighbors in a list
def pairwise(xs):
    return [(xs[idx],xs[idx+1]) for idx in range(len(xs)-1)]

# Return list of pairs of pairwise neighbor INDEXES in a list
def pairwise_indexes(xs):
    return [(idx,idx+1) for idx in range(len(xs)-1)]

# Merge inner lists by initial index (must supply >= 2 indexes)
def mergeInnerLists(xss, idxs):
    assert len(idxs) >= 2
    mergeIntoIdx = min(idxs)
    seenMergeIntoIdx = False
    res = []
    for (jdx,xs) in enumerate(xss):
        if jdx in idxs:
            if seenMergeIntoIdx:
                res[mergeIntoIdx] += xs
            else:
                res += [xs]
                seenMergeIntoIdx = True
        else:
            res += [xs]
    return res

# Explode a list into a list of singleton lists
def explode(xs):
    return [[x] for x in xs]
