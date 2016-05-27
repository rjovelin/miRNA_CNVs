# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 14:22:37 2015

@author: RJovelin
"""


# use this script to split the targetscan sequence input file into several input files

import os
import sys


# usage python3 split_targetscan_input.py species directory Ngenes

# get species from command line
species = sys.argv[1]

# get destination directory from command
directory = sys.argv[2]

# get the number of genes to add in each inputfile
Ngenes = int(sys.argv[3])

# get targetscan input file
targetscan_input = species + '_5UTR_valid_chromos_targetscan.txt'

# create a dict {gene: list} from inputfile
gene_seq = {}

# create a list with all genes in same order as in the input file
genes = []

# open targetscan input file for reading
infile = open(targetscan_input, 'r')
# loop over file
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split()
        # populate dict
        gene_seq[line[0]] = '\t'.join(line)
        # add gene to list
        genes.append(line[0])
# close file after reading
infile.close()

# loop over list of genes
for i in range(0, len(genes) + Ngenes, Ngenes):
    print(i)
    # open targetscan inputfile
    newfile = open(directory + species + '_5UTR_valid_chromos_targetscan_' + str(i) + '.txt', 'w')
    # loop over genes within this slice of the list
    for gene in genes[i:i + Ngenes]:
        newfile.write(gene_seq[gene] + '\n')
    # close file after writing
    newfile.close()