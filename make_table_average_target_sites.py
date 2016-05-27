# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 14:48:14 2015

@author: RJovelin
"""


# usage python3 make_table_average_target_sites.py [3UTR/5UTR/CDS] [targetscan/miranda] [True/False] [long_CNVs/all_CNVs]

# write to file average number of targets, sequence length for CNV and non-CNV genes per species

from CNV_miRNAs import *
import os
import sys
import numpy as np
from scipy import stats
import math

# get the region to consider to predict target sites [3UTR or 5UTr or CDS]
domain = sys.argv[1]
print(domain)

# get the predictor [targetscan or miranda]
predictor = sys.argv[2]
print(predictor) 

# get the option to keep genes on all chromos (False) or only on assembled 
# nuclear chromosomes only from the command
keep_valid_chromos = sys.argv[3]
if keep_valid_chromos == 'True':
    keep_valid_chromos = True
    chromos = 'valid_chromos'
elif keep_valid_chromos == 'False':
    keep_valid_chromos = False
    chromos = 'all_chromos'
print(keep_valid_chromos, chromos)

# check if all chromos (including unplaced, unlocated, and MT) are used
# or if only valid chromos are used 
# note: only CNVs on valid chromos are reported in DGV, so if all chromos are
# used it may introduce a bias by calling non CNV genes genes that cannot be detected

# get the type of CNVs to consider from the command
long_CNV = sys.argv[4]
if long_CNV == 'all_CNVs':
    cnv_length = 'CNV_all_length'
elif long_CNV == 'long_CNVs':
    cnv_length = 'CNV_greater_1Kb'
print(long_CNV, cnv_length)

# get outputfile
outputfile = 'Average_targets_' + domain + '_' + chromos + '_' + cnv_length + '_' + predictor + '.txt'
print(outputfile)

# open file for writing
newfile = open(outputfile, 'w')

# write header to file
newfile.write('\t'.join(['Species', 'N_CNV_genes', 'CNV_mean_targets', 'CNV_SEM_targets',  
                        'N_nonCNV_genes', 'nonCNV_mean_targets', 'nonCNV_SEM_targets',
                        'P_diff_targets', 'CNV_mean_normalized_targets', 'CNV_SEM_normalized_targets',
                        'nonCNV_mean_normalized_targets', 'nonCNV_SEM_normalized_targets', 'P_diff_normalized_targets',
                        'CNV_mean_seq_length', 'CNV_SEM_seq_length', 'nonCNV_mean_seq_length',
                        'nonCNV_SEM_seq_length', 'P_seq_length', 'Spearman_rho_targets_X_length', 'P_targets_X_length']) + '\n')                        


# create a lamda function to transform value into string
Gstr = lambda x: str(x)

# make a dictionary of species names : species code
species_names = {'H_sapiens': 'Hsa',  'P_troglodytes': 'Ptr', 'M_mulatta': 'Mmul',
                 'M_musculus': 'Mmus', 'R_norvegicus': 'Rno', 'B_taurus': 'Bta',
                 'C_familiaris': 'Cfa', 'G_gallus': 'Gga'}

# loop over species names
for species in species_names:
    print(species)
    
    if species == 'H_sapiens':
        # get summary table for DGV 2015 release
        infile = 'H_sapiens_' + domain + '_summary_' + predictor + '_' + chromos + '_' + cnv_length + '_GRCh37_2015.txt'
    else:
        # get summary table 
        infile = species + '_' + domain + '_summary_' + predictor + '_' + chromos + '_' + cnv_length + '.txt'
    print(infile)
    
    # parse the summary table into a list
    regulation = compare_miRNA_regulation(infile)
    
    # write regulation to file
    newfile.write(species + '\t')
    newfile.write('\t'.join(list(map(Gstr, regulation))) + '\n')
    
    print('done writing regulation for {0}'.format(species))

# close file after writing
newfile.close()