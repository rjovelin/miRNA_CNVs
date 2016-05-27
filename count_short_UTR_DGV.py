# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 10:24:09 2015

@author: RJovelin
"""


# usage count_short_3UTR_DGV.py [True/False] [long_CNVs/all_CNVs]


# determine the number of genes with short (< 7b) 3'UTR in each species 

from CNV_miRNAs import *
import sys

# get the option to keep genes on all chromos (False) or only on assembled 
# nuclear chromosomes only (True) from the command

# check if all chromos (including unplaced, unlocated, and MT) are used
# or if only valid chromos are used 
# note: only CNVs on valid chromos are reported in DGV, so if all chromos are
# used it may introduce a bias by calling non CNV genes genes that cannot be detected

keep_valid_chromos = sys.argv[1]
if keep_valid_chromos == 'True':
    keep_valid_chromos = True
    chromos = 'valid_chromos'
elif keep_valid_chromos == 'False':
    keep_valid_chromos = False
    chromos = 'all_chromos'
print(keep_valid_chromos, chromos)   


# get the option to call CNV based on CNV length from the command
long_CNV = sys.argv[2]
if long_CNV == 'all_CNVs':
    cnv_length = 'CNV_all_length'
elif long_CNV == 'long_CNVs':
    cnv_length = 'CNV_greater_1Kb'

# get UTR file
UTR_file = 'H_sapiens_3UTR_length_' + chromos + '.txt'
print(UTR_file)

# make a list of CNV files
CNV_files = ['H_sapiens_GRCh37_2013-05_' + cnv_length + '_' + chromos + '.txt',
             'H_sapiens_GRCh37_2013-07_' + cnv_length + '_' + chromos + '.txt',
             'H_sapiens_GRCh37_2014_' + cnv_length + '_' + chromos + '.txt',
             'H_sapiens_GRCh37_2015_' + cnv_length + '_' + chromos + '.txt']
             
# get outputfile
outputfile =  'Gene_Counts_DGV_release_short_3UTR_' + cnv_length + '_' + chromos + '.txt'             
print(outputfile)             
             
# open file for writing
newfile = open(outputfile, 'w')             
             
# write headers to file
newfile.write('# Number of genes with short (< 7 bp) 3\'UTR\n')
newfile.write('\t'.join(['DGV_release', 'Total', 'CNV', 'non-CNV']) + '\n')

# loop over CNV files
for filename in CNV_files:
    print(filename)
    # sort genes based on CNV status
    CNV_status = sort_genes_CNV_status(filename)
    print(len(CNV_status))
    # sort genes based on 3' UTR length
    UTR_length = sort_genes_3UTR_length(UTR_file)
    print(len(UTR_length))
    # get release version
    release_version = filename[filename.index('GRCh37'): filename.index('_CNV')]
    print(release_version)
    # count total number of genes with short 3'UTR
    total_short = 0
    # count CNV genes with short UTR
    cnv_short = 0
    # count non-CNV genes with short UTR
    non_cnv_short = 0    
    # loop over genes in UTR_length
    for gene in UTR_length:
        if UTR_length[gene] == 'short':
            total_short += 1
            # check CNV status
            if CNV_status[gene] == 'CNV':
                cnv_short += 1
            elif CNV_status[gene] == 'not_CNV':
                non_cnv_short += 1
            
    print('total short', total_short)
    print('CNV short', cnv_short)
    print('non_CNV', non_cnv_short)
    
    assert total_short == cnv_short + non_cnv_short, 'sum cnv and non-cnv short is not equal to total short'
        
    # write results to file
    newfile.write('\t'.join([release_version, str(total_short), str(cnv_short), str(non_cnv_short)]) + '\n')
    
# close file after writing
newfile.close()

