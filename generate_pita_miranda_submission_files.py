# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 00:00:35 2015

@author: Richard
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
files = [i for i in os.listdir() if species in i and domain in i and 'miranda' in i]

# loop over input files
for filename in files:
    # get number from filename
    num = filename[filename.rindex('_') + 1: filename.index('.txt')]
    # get outputfile
    outputfile = species + '_run_miranda_' + domain + '_' + num + '.sh'
    # open submission file for writing
    newfile = open(outputfile, 'w')
    # predict target sites with miranda
    # get command
    command = '../miranda ' + '../' + species + '_mature.txt ' + \
    filename + ' -en -10 -sc 140 -strict -out ' + \
    species + '_' + domain + '_valid_chromos_predicted_miranda_' + num + '.txt'
    # write command to file
    newfile.write(command)
    # close file after writing
    newfile.close()
    
