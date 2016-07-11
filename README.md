# Clustal-MUSCLE-Iterator

A script to iterate through the parameters of Clustal/MUSCLE to generate the best sum of pair score.

This script will generate mega alignment configurations (.mao) and perform Clustal or MUSCLE alignments using the
megasoftware toolkit.

## Getting Started

### Prerequisites

[Python3](https://www.python.org/downloads)

[BioPython](http://biopython.org/wiki/Download)
- Used for parsing fasta formatted files and supplies substituion matrices

[Numpy](https://www.scipy.org/install.html)
- Dependency for BioPython

[MEGA](http://www.megasoftware.net/)
- Download the MegaCC (command line) version for your os
- Used to generate multiple sequence alignments

### Setup
- Install Python3, BioPython, Numpy
- Clone this repo
- Download MEGA and extract it to the mega folder (or remember the path)
