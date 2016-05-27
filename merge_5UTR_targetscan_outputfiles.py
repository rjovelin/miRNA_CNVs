# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 12:22:00 2015

@author: RJovelin
"""


# use this script to merge the targetscan outputfiles for 5UTRs

import os
import sys

# usage


# get the outputfile from the command
outputfile = sys.argv[1]


# make a list of targetscan output files
targetscan_out = [i for i in os.listdir() if '5UTR' in i and 'predicted' in i]

# create a list with num in file and corresponding filename
nums = {}
for filename in targetscan_out:
    number = int(filename[filename.rindex('_') + 1: filename.index('.txt')])
    nums[number] = filename
# make a list of numbers in file
numbers = [i for i in nums]
# sort list
numbers.sort()

# open outputfile for writing
newfile = open(outputfile, 'w')

# write header
header = ['a_Gene_ID', 'miRNA_family_ID species_ID', 'MSA_start', 'MSA_end UTR_start', 'UTR_end Group_num',
          'Site_type', 'miRNA in this species', 'Group_type', 'Species_in_this_group', 'Species_in_this_group_with_this_site_type']
newfile.write('\t'.join(header) + '\n')


# loop over sorted numbers, grab corresponding file
for i in numbers:
    filename = nums[i]
    # open file for reading
    infile = open(filename, 'r')
    # loop over file
    for line in infile:
        # skip header and empty lines
        if not line.startswith('a_Gene_ID') and line.rstrip() != '':
            # write line to newfile
            newfile.write(line)
    # close infile after reading
    infile.close()
    
# close file after writing
newfile.close()
        


