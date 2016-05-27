# -*- coding: utf-8 -*-
"""
Created on Sun Oct  4 12:39:04 2015

@author: Richard
"""

# usage python3 make_valid_chromo_files.py [directory]

import sys
import os

directory = sys.argv[1]


# save each chromosome of the nuclear genome in a file

# make a list of fasta files for each chromosome
files = [i for i in os.listdir(directory) if 'chr' in i and 'chrMT' not in i]

# make a set of valid chromosomes
chromos = set()

# loop over files
for filename in files:
    # open file for reading
    infile = open(directory + filename, 'r')
    # get header
    header = infile.readline().rstrip().split('|')
    chromos.add(header[3])
    # close file
    infile.close()
    # open file for writing
    newfile = open(directory[:directory.index('/')] + '_valid_chromos.txt', 'w')
    # loop over set of chromos
    for LG in chromos:
        newfile.write(LG + '\n')
    # close file
    newfile.close()
    

    