# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 14:21:28 2016

@author: RJovelin
"""

# use this script to make a file of chromosome names matched to scaffold/contig names

import os



# make a list of fasta files
files = [i for i in os.listdir() if i[-3:] == '.fa']

# make a dictionary of chromosome: scaffold
chromos = {}

# loop over fasta files
for filename in files:
    infile = open(filename)
    # look for sequence headers
    for line in infile:
        if line.startswith('>'):
            line = line.rstrip().split('|')
            scaffold = line[3]
            LG = line[-1]
            LG = LG[:LG.index(',')]
            LG = LG.replace('Homo sapiens', '')
            LG = LG.replace(' ', "")
            LG = LG.replace('omosome', '')
            if LG in chromos:
                chromos[LG].append(scaffold)
            else:
                chromos[LG] = [scaffold]
    infile.close()

                
                
