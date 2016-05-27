# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 20:16:51 2015

@author: Richard
"""


# usage python3 compare_DVG_release_valid_chromos.py [3UTR/5UTR/CDS] [targetscan/miranda] [all_CNVs/long_CNVs]

# write to file average number of targets, sequence length for CNV and non-CNV genes for each DVG release


from CNV_miRNAs import *
import os
import sys
import numpy as np
from scipy import stats
import math

# get the region to consider to predict target sites [3UTR or 5UTr or CDS]
domain = sys.argv[1]
print(domain)

# get program used to predict target sites [targetscan or miranda]
predictor = sys.argv[2]
print(predictor)

# get the option to call a CNV if CNV length > 1 Kb (long_CNVs)
# or to include all CNVs regardless of length (all_CNVs)
long_CNV = sys.argv[3]
print(long_CNV)

# check if all CNVs are considered or only CNVs > 1 Kb
if long_CNV == 'all_CNVs':
    cnv_length = 'CNV_all_length'
elif long_CNV == 'long_CNVs':
    cnv_length = 'CNV_greater_1Kb'

# make a list of CNV gene files
summary_files = ['H_sapiens_' + domain + '_summary_' + predictor + '_valid_chromos_' + cnv_length + '_GRCh37_2013-05.txt',
                 'H_sapiens_' + domain + '_summary_' + predictor + '_valid_chromos_' + cnv_length + '_GRCh37_2013-07.txt',
                 'H_sapiens_' + domain + '_summary_' + predictor + '_valid_chromos_' + cnv_length + '_GRCh37_2014.txt',
                 'H_sapiens_' + domain + '_summary_' + predictor + '_valid_chromos_' + cnv_length + '_GRCh37_2015.txt']

# get outputfile
outputfile = 'DGV_release_average_targets_' + predictor + '_' + domain + '_valid_chromos_' + cnv_length + '.txt'
print(outputfile)

# open file for writing
newfile = open(outputfile, 'w')

# write header to file
newfile.write('\t'.join(['Release', 'N_CNV_genes', 'CNV_mean_targets', 'CNV_SEM_targets',  
                        'N_nonCNV_genes', 'nonCNV_mean_targets', 'nonCNV_SEM_targets',
                        'P_diff_targets', 'CNV_mean_normalized_targets', 'CNV_SEM_normalized_targets',
                        'nonCNV_mean_normalized_targets', 'nonCNV_SEM_normalized_targets', 'P_diff_normalized_targets',
                        'CNV_mean_seq_length', 'CNV_SEM_seq_length', 'nonCNV_mean_seq_length',
                        'nonCNV_SEM_seq_length', 'P_seq_length', 'Spearman_rho_targets_X_length', 'P_targets_X_length']) + '\n')                        

# create a lamda function to transform value into string
Gstr = lambda x: str(x)


# loop over summary table, compare CNV and non-CNV genes
for filename in summary_files:
    print(filename)
    
    # parse the summary table into a list
    regulation = compare_miRNA_regulation(filename)
    
    # get release number
    release_version = filename[filename.index('GRCh37') : filename.index('.txt')]
    print(release_version)
        
    # write regulation to file
    newfile.write(release_version + '\t')
    newfile.write('\t'.join(list(map(Gstr, regulation))) + '\n')
    
    print('done writing regulation for {0}'.format(release_version))

# close file after writing
newfile.close()