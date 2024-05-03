# Localized Enumerator

Code to accompany the paper "Geometric enumeration of localized DNA strand displacement reaction networks".

## Description

- This package implements a prototype localized variant of the previously reported geometric reaction enumeration framework for domain level DNA strand displacement reactions.
- The code demonstrates that the geometric framework described in that paper can be used to enumerate localized reactions and implements the examples presented in that paper, using a structure sampling approach for constraint solving.
- For further details please refer to the paper.

## Getting Started
- The package contains this README file, a LICENSE file and the src folder.
- The src folder contains the code for reaction enumeration.

### Dependencies

* Python 3.x
* ply
* Jupyter Notebook (optional for displaying the graphs)
* Graphviz library (optional for displaying the graphs)

### Executing programs

- The jupyter notebooks in the src directory contain code for running the reactions illustrated in the manuscript.
- Those examples, and the accompanying paper, contain further details of the input syntax, which is briefly outlined below.
- Some of the classes define graphical representations that can be used for visual debugging and development.
  The graph visualization relies on the "graphviz" library and associated command-line tool being installed.

### Input/output formats

- The input DNA strands are provided using a variant of the previously-reported process calculus syntax, which includes representations of tile species and tether locations. The lengths of the domains are provided separately. An example: 
```
# Domain length information:
domainLengthStr = 'toeholdDomain spcr1 length 5 toeholdDomain spcr2 length 5 longDomain s length 12 toeholdDomain a0 length 6 toeholdDomain f length 6 toeholdDomain x length 6 longDomain y length 12'

# Details of DNA strands through process calculus syntax
s = '( <s a0^> | [[<tether(0,0) spcr1 a0^* s*!i1 f^ s!i1> | <tether(10.88,0) spcr2 x^* s*!i3 y^ s!i3> ]] | <s!i2 x^ s*!i2 f^*> | <s*!i4 y^*> | <s!i4> )'
```
- Output can be produced in the following two formats:
  - Command line:
    - Lists all the possible reactions that can occur and
    - The strand graph information details of the species that are part of the reaction.
  - Jupyter notebook: 
    - Lists all the possible reactions and 
    - Displays the graphical version of species and can be saved as Jupyter notebooks in a html or pdf format.

## Authors

Matthew R. Lakin
Sarika Kumar

## Acknowledgments

This material is based upon work supported by the National Science Foundation under Grants 1518861 and 1814906.
