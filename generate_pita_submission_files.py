# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 14:16:06 2015

@author: RJovelin
"""

# usage generate_pita_miranda_submission_files.py species domain [3UTR/5UTR/CDS] 

# use this script to generate submission file to run targetscan for a single file
# run this script in the same directory with input sequence files


import os
import sys

# get species from command line
species = sys.argv[1]
print(species)

domain = sys.argv[2]
print(domain)

# make a list of sequence input files
files = [i for i in os.listdir() if species in i and 'mature' in i]

# loop over input files
for filename in files:
    # get the mirna name from the file name
    mirna = filename[filename.rindex('_'): -4]
    # replace empty space
    if ' ' in mirna:
        mirna = mirna.replace(' ', '-')
    # get outputfile
    outputfile = species + '_run_pita_' + domain + '_' + mirna + '.sh'
    # open submission file for writing
    newfile = open(outputfile, 'w')
    # predict target sites with PITA
    # get command
    command = 'perl pita_prediction.pl -mir ' + filename + \
    ' -utr ' + '../' + species + '_' + domain + '_valid_chromos_fasta.txt' + \
    ' -l 6-8 -prefix ' + species + '_' + domain + '_valid_chromos_predicted_pita_' + mirna + '.txt'
    # write command to file
    newfile.write(command)
    # close file after writing
    newfile.close()
    
