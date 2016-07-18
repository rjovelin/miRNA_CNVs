# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 14:42:11 2015

@author: RJovelin
"""

# usage 


# usage: python3 make_human_CNV_file.py [long/all] 

import os
import sys
from CNV_miRNAs import *


# get CNV_size from the command
CNV_size = sys.argv[1]

# check CNV_size
if CNV_size == 'all':
    cnv_length = 'all_length'
elif CNV_size == 'long':
    cnv_length = 'greater_1Kb'


DGV_files = ['GRCh37_hg19_variants_2015-07-23.txt', 'GRCh37_hg19_variants_2014-10-16.txt',
             'GRCh37_hg19_variants_2013-07-23.txt', 'GRCh37_hg19_variants_2013-05-31.txt']  

# loop over DGV files
for filename in DGV_files:
    print(filename)
    # make a set of CNV genes
    CNV_genes = get_human_CNV_genes(filename, CNV_size)
    print('H_sapiens CNV genes', len(CNV_genes))
    
    # get outputfile
    if '2013' in filename:
        outputfile = 'human_CNV_' + cnv_length + \
        '_GRCh37_' + filename[filename.rindex('_') + 1: -7] + '.txt'
    else:
        outputfile = 'human_CNV_' + cnv_length + \
        '_GRCh37_' + filename[filename.rindex('_') + 1: -10] + '.txt'
    print(outputfile)

    # open file for writing
    newfile = open(outputfile, 'w')
    
    # loop over CNV genes, write gene to file
    for gene in CNV_genes:
        newfile.write(gene + '\n')
    # close file after writing
    newfile.close()