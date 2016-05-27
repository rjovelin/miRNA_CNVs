# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

# usage run_targetscan_prediction_5UTR.py walltime


# use this script to submit jobs to predict target sites in 5'UTR 
# with input sequence file split into several files

import os
import sys

# get walltime from the command
walltime = sys.argv[1]

# make alist of submission files
submission_files = [i for i in os.listdir() if '5UTR' in i and '.sh' in i]

# loop over file, submit file to cluster
for filename in submission_files:
    os.system('qsub -d `pwd` -l walltime=' + str(walltime) + ':0:0,nodes=1:ppn=1 ' + filename)
    
    
