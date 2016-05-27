# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 14:37:06 2015

@author: RJovelin
"""

# usage generate_targetscan_submission_files.py species

# use this script to generate submission file to run targetscan for a single file

import os
import sys

# get species from command line
species = sys.argv[1]

# make a list of targetscan sequence input files
targetscan_files = [i for i in os.listdir() if species in i and '5UTR' in i and 'targetscan' in i]

# loop over input files
for filename in targetscan_files:
    # get number from filename
    num = filename[filename.rindex('_') + 1: filename.index('.txt')]
    # open submission file for writing
    newfile = open(species + '_run_targetscan_5UTR_' + num + '.sh', 'w')
    newfile.write('perl targetscan_60.pl ' + species + '_miRFam_targetscan.txt ' + filename + ' ' + species + '_5UTR_valid_chromos_predicted_' + num + '.txt')
    newfile.close()
