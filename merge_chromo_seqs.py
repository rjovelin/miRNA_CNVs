# -*- coding: utf-8 -*-
"""
Created on Sun Sep 27 09:12:59 2015

@author: Richard
"""

import os
from CNV_miRNAs import *

# make a list of directories
directories = ['B_taurus/', 'C_familiaris/', 'G_gallus/',  'H_sapiens/',
               'M_mulatta/',  'M_musculus/', 'P_troglodytes/', 'R_norvegicus/']
               
# loop over directories
for i in range(len(directories)):
    print(directories[i])
    # create a dictionary to hold the sequences
    genome = {}    
    # loop over files in directory
    for filename in os.listdir(directories[i]):
        print(filename)
        # check if file contains unplaced or unlocalized sequences
        if 'unplaced' in filename or 'unlocalized' in filename:
            print('use convert_fasta')
            # multiple sequences in file, use convert_fasta
            fasta_seq = convert_fasta(directories[i] + filename)
        else:
            print('use read_genome')
            # single sequence, use read_chromo
            fasta_seq = read_genome(directories[i] + filename)
        # populate genome dict
        for chromo in fasta_seq:
            genome[chromo] = fasta_seq[chromo]
    # open file for writing
    newfile = open(directories[i][:-1] + '_genome.txt', 'w')
    print(directories[i][:-1] + '_genome.txt')
    # loop over chromo in genome, write sequence in fasta to file
    for chromo in genome:
        newfile.write('>' + chromo + '\n')
        newfile.write(genome[chromo] + '\n')
    # close file
    newfile.close()
  

