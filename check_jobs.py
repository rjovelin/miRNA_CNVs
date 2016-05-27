# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 21:38:26 2015

@author: Richard
"""


# usage check_jobs.py


import os

# use this script to check if jobs were killed

files = [i for i in os.listdir() if '.sh.e' in i]

# loop over error files
for filename in files:
    # open file for reading
    infile = open(filename, 'r')
    # get content
    content = infile.read()
    # close file
    infile.close()
    # check if job was killes
    if 'job killed' in content:
        print(filename, content, 'jobs killed prematurely', sep = '\n')
        break
    

        

