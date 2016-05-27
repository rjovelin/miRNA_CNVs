# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 16:48:14 2015

@author: RJovelin
"""

# convert targetscan input sequence into a fasta file for input by miRanda and PITA

# usage convert_targetscan_seq_input_to_fasta.py [3UTR/5UTR/CDS] [True/False]


import sys
import os

# get domain from command
domain = sys.argv[1]
print(domain)

# get the option to keep assembled chromos or all chromos
keep_valid_chromos = sys.argv[2]
if keep_valid_chromos == 'True':
    keep_valid_chromos = True
    chromos = 'valid_chromos'
elif keep_valid_chromos == 'False':
    keep_valid_chromos = False
    chromos = 'all_chromos'
print(keep_valid_chromos, chromos)
    
# make a list of species
species_names = ['H_sapiens', 'P_troglodytes', 'M_mulatta', 'M_musculus',
                 'R_norvegicus', 'B_taurus', 'C_familiaris', 'G_gallus']

# loop over species
for species in species_names:
    print(species)
    # get targetscan input file
    targetscan_input = species + '_' + domain + '_' + chromos + '_targetscan.txt'
    print(targetscan_input)
    # create a dict of {gene : sequence} pairs
    gene_seq = {}
    # open file for reading
    infile = open(targetscan_input, 'r')
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            gene_seq[line[0]] = line[2]
    infile.close()
    # open file for writing
    outputfile = species + '_' + domain + '_' + chromos + '_fasta.txt'
    newfile = open(outputfile, 'w')
    # loop over genes
    for gene in gene_seq:
        newfile.write('>' + gene + '\n')
        newfile.write(gene_seq[gene] + '\n')
    newfile.close()

 
 

