# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 13:17:29 2015

@author: RJovelin
"""


# usage run_miranda_prediction_5UTR.py walltime domain [3UTR/5UTR/CDS] 

# run this script in the directory where submission files are located

# use this script to submit jobs to predict target sites by miranda 
# with input sequence file split into several files

import os
import sys

# get walltime from the command
walltime = sys.argv[1]

# get domain
domain = sys.argv[2]

# make alist of submission files
submission_files = [i for i in os.listdir() if domain in i and 'miranda' in i and '.sh' in i]

# loop over file, submit file to cluster
for filename in submission_files:
    os.system('qsub -d `pwd` -l walltime=' + str(walltime) + ':0:0,nodes=1:ppn=1 ' + filename)
    
    
