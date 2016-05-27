# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 13:57:44 2015

@author: RJovelin
"""


# use this script to split the mature miRNAs with a single miRNA per file

import os
import sys
from CNV_miRNAs import *

# usage python3 split_mature_input.py species destination [destination directory] 

# get species from command line
species = sys.argv[1]

# get destination directory from command
destination = sys.argv[2]

# get mature fasta file
mature_fasta = convert_fasta(species + '_mature.txt')

# loop over mature mirnas
for mirna in mature_fasta:
    # replace empty space
    if ' ' in mirna:
        mirname = mirna.replace(' ', '-')
    # open mature input file
    newfile = open(destination + species + '_mature_' + mirname + '.txt', 'w')
    # write mirna and mirna to file in fasta
    newfile.write('>' + mirna + '\n')
    newfile.write(mature_fasta[mirna] + '\n')
    newfile.close()
    
    
