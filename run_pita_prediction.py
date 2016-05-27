# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 14:31:05 2015

@author: RJovelin
"""

# usage run_pita_prediction_5UTR.py walltime domain [3UTR/5UTR/CDS] directory

# run this script in the directory where submission files are located

# use this script to submit jobs to predict target sites by pita 
# with mature input sequence file split into several files

import os
import sys

# get walltime from the command
walltime = sys.argv[1]

# get domain
domain = sys.argv[2]

# get directory where results sbmission files and results are stored 
directory = sys.argv[3]

# make a list of submission files that have already been run
already_launched = [i for i in os.listdir(directory) if i[-3:] == '.sh']

# make alist of submission files
submission_files = [i for i in os.listdir() if domain in i and 'pita' in i and '.sh' in i]

# run only 500 jobs at a time (maximum number of jobs)
j = 400
# loop over file, submit file to cluster
for filename in submission_files:
    print(j)
    # check if maximum jobs is reached
    if j != 0:
        # check if job has previously been launched
        if filename not in already_launched:
            # run job
            os.system('qsub -d `pwd` -l walltime=' + str(walltime) + ':0:0,nodes=1:ppn=1 ' + filename)
            # update counter
            j -= 1
    
    
