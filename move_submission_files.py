# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 10:11:39 2015

@author: RJovelin
"""

# use this script to move submission files and corresponding err files
# to a directory once the jobs are done

# usage move_submission_files.py directory

import os
import sys

# get directory from command

directory  = sys.argv[1]

# make a list of error files
errfiles = [i for i in os.listdir() if '.sh.e' in i]

for filename in errfiles:
    # get corresponding submission file
    subfile = filename[:filename.index('.sh.e')] + '.sh'
    # move submission file to directory
    os.system('mv ' + subfile + ' ' + directory)
    # move error file to directory
    os.system('mv ' + filename + ' ' + directory)