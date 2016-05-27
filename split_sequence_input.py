# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 23:05:52 2015

@author: Richard
"""


# use this script to split the targetscan sequence input file into several input files

import os
import sys
from CNV_miRNAs import *

# usage python3 split_sequence_input.py species destination [destination directory] domain [CDS/3UTR/5UTR] predictor [miranda/pita] Ngenes

# get species from command line
species = sys.argv[1]

# get destination directory from command
destination = sys.argv[2]

# get domain
domain = sys.argv[3]

# get predictor
predictor = sys.argv[4]

# get the number of genes to add in each inputfile
Ngenes = int(sys.argv[5])

# get targetscan input file
fasta_seq = species + '_' + domain + '_valid_chromos_fasta.txt'

# create a dict {gene: sequence} from inputfile
gene_seq = convert_fasta(fasta_seq)

# create a list with all genes in same order as in the input file
genes = []

# open seq input file for reading
infile = open(fasta_seq, 'r')
# loop over file
for line in infile:
    if line.startswith('>'):
        line = line.rstrip()
        genes.append(line[1:])
# close file after reading
infile.close()

# loop over list of genes
for i in range(0, len(genes) + Ngenes, Ngenes):
    print(i)
    # open targetscan inputfile
    newfile = open(destination + species + '_' + domain + '_valid_chromos_' + predictor + '_' + str(i) + '.txt', 'w')
    # loop over genes within this slice of the list
    for gene in genes[i:i + Ngenes]:
        newfile.write('>' + gene + '\n')
        newfile.write(gene_seq[gene] + '\n')
    # close file after writing
    newfile.close()