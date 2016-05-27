# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 15:48:36 2015

@author: RJovelin
"""


# use this script to merge the miranda outputfiles from split job results

import os
import sys

# usage python3 merge_miranda_predictions.py species domain [3UTR/5UTR/CDS] outputfile

# run this script in directory with miranda outputfiles

# get species from command
species = sys.argv[1]

# get sequence domain
domain = sys.argv[2]

# get the outputfile from the command
outputfile = sys.argv[3]


# make a list of targetscan output files
miranda_out = [i for i in os.listdir() if domain in i and 'predicted' in i]

# create a list with num in file and corresponding filename
nums = {}
for filename in miranda_out:
    number = int(filename[filename.rindex('_') + 1: filename.index('.txt')])
    nums[number] = filename
# make a list of numbers in file
numbers = [i for i in nums]
# sort list
numbers.sort()

# open outputfile for writing
newfile = open(outputfile, 'w')

# loop over sorted numbers, grab corresponding file
for i in numbers:
    filename = nums[i]
    # open file for reading
    infile = open(filename, 'r')
    # loop over file
    for line in infile:
        # record only lines corresponding to complete scan
        # of a given mirna for a given gene
        if line.startswith('>>'):
            # write line to merged file
            newfile.write(line)
    # close infile after reading
    infile.close()
    
# close file after writing
newfile.close()
        


