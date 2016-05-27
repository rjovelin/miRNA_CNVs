# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 11:36:03 2015

@author: RJovelin
"""


# usage check_jobs_2.py


import os

# use this script to check if jobs were killed

files = [i for i in os.listdir() if '.sh.e' in i]

# loop over error files
for filename in files:
    # open file for reading
    infile = open(filename, 'r')
    # check first line
    line = infile.readline().rstrip()
    if line.startswith('/RQusagers/'):
        # get next line
        line = infile.readline().rstrip()
        if line.startswith('/RQusagers/'):
            line = infile.readline().rstrip()
    
    # check if line is empty or has error message
    if line != '':
        print(filename, line, sep = '\n')
        # close file
        infile.close()
        break
    

        

