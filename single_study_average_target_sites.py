# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 11:08:48 2015

@author: RJovelin
"""


# usage single_study_average_target_sites.py [parameters]
# [3UTR/5UTR/CDS] choose the region to analyse
# [targetscan/miranda] choose the target predictor
# [True/False] use valid chromos (True) or all chromos
# [long_CNVs/all_CNVs] use all CNVs or CNVs > 1 Kb
# CNV_file Choose DGV CNV file
# minimum_cnv minimum number of cnv genes in single study


from CNV_miRNAs import *
import os
import sys
import numpy as np
from scipy import stats
import math

# get the region to consider to predict target sites [3UTR or 5UTr or CDS]
domain = sys.argv[1]
print(domain)

# get predictor from command [targetscan or miranda]
predictor = sys.argv[2]
print(predictor)

# get the option to keep genes on all chromos (False) or only on assembled 
# nuclear chromosomes only (True) from the command
keep_valid_chromos = sys.argv[3]
if keep_valid_chromos == 'True':
    keep_valid_chromos = True
    chromos = 'valid_chromos'
elif keep_valid_chromos == 'False':
    keep_valid_chromos = False
    chromos = 'all_chromos'
print(keep_valid_chromos, chromos)


# get the option to call a CNV if CNV length > 1 Kb (long_CNVs)
# or to include all CNVs regardless of length (all_CNVs)
long_CNV = sys.argv[4]
if long_CNV == 'all_CNVs':
    cnv_length = 'CNV_all_length'
    CNV_size = 'all'
elif long_CNV == 'long_CNVs':
    cnv_length = 'CNV_greater_1Kb'
    CNV_size = 'long'
print(long_CNV, cnv_length, CNV_size)

# get the CNV file from the command line
# ie. all_CNV or CNVs > 1 Kb can be used for any DGV release
CNV_file = sys.argv[5]
print(CNV_file)

# get the minimum number of cnv genes that a single study must have
minimum_cnv = int(sys.argv[6])
print(minimum_cnv)


# check if all chromos (including unplaced, unlocated, and MT) are used
# or if only valid chromos are used 
# note: only CNVs on valid chromos are reported in DGV, so if all chromos are
# used it may introduce a bias by calling non CNV genes genes that cannot be detected

# get the file with 3'UTR length
UTR_file = 'H_sapiens_3UTR_length_' + chromos + '.txt'
print(UTR_file)
# get targetscan sequence input file
targetscan_seq_input_file = 'H_sapiens_' + domain + '_' + chromos + '_targetscan.txt'
print(targetscan_seq_input_file)

# get the outputfile with predicted target sites
predicted_target_file = 'H_sapiens_' + domain + '_' + chromos + '_predicted_sites_' + predictor + '.txt'    
print(predicted_target_file)
# make a dictionary with {gene :[targets, seq_length, normalized_targets]}

# check predictor
if predictor == 'targetscan':
    predicted_targets = parse_targetscan_output(targetscan_seq_input_file, predicted_target_file, 'all')
elif predictor == 'miranda':
    predicted_targets = parse_miranda_output(targetscan_seq_input_file, predicted_target_file, 'all')

# get release version
if '2013' in CNV_file:
    release_version = 'GRCh37_' + CNV_file[CNV_file.rindex('_') + 1: -7]
else:
    release_version = 'GRCh37_' + CNV_file[CNV_file.rindex('_') + 1: -10]
print(release_version)

# get outputfile
outputfile = 'H_sapiens_single_study_targets_' + domain + '_' + cnv_length + '_' + chromos + '_' + predictor + '_' + release_version + '.txt'
print(outputfile)

# open file for writing
newfile = open(outputfile, 'w')

# write header to file
newfile.write('\t'.join(['Study', 'N_CNV_genes', 'CNV_mean_targets', 'CNV_SEM_targets',  
                        'N_nonCNV_genes', 'nonCNV_mean_targets', 'nonCNV_SEM_targets',
                        'P_diff_targets', 'CNV_mean_normalized_targets', 'CNV_SEM_normalized_targets',
                        'nonCNV_mean_normalized_targets', 'nonCNV_SEM_normalized_targets', 'P_diff_normalized_targets',
                        'CNV_mean_seq_length', 'CNV_SEM_seq_length', 'nonCNV_mean_seq_length',
                        'nonCNV_SEM_seq_length', 'P_seq_length', 'Spearman_rho_targets_X_length', 'P_targets_X_length']) + '\n')                        

# get the dictionaries of {reference: pubmedid} for studies reported in the CGV
references = get_DGV_references(CNV_file)

# create a lamda function to transform value into string
Gstr = lambda x: str(x)

# loop over study
for study in references:
    print(study)
    # get the set of CNV genes corresponding to that study
    CNV_genes = get_human_CNV_genes_single_study(CNV_file, study, CNV_size)
    print('# CNV genes', len(CNV_genes)) 
    
    # check that study includes minimum number of cnv genes
    if len(CNV_genes) >= minimum_cnv:
        # make temporary cnv_file with CNV genes extracted from DGV
        tempfile = open('Temp_cnv_file.txt', 'w')
        # dump all CNV genes
        for gene in CNV_genes:
            tempfile.write(gene + '\n')
        tempfile.close()
    
        # get CNV gene status
        CNV_status = get_genes_CNV_status('H_sapiens.gff3', 'H_sapiens_genome.txt', 'H_sapiens_valid_chromos.txt', keep_valid_chromos, 'Temp_cnv_file.txt')    
        print('genes with CNV status', len(CNV_status))
    
        # make a temporary file with CNV status of all genes     
        tempfile = open('Temp_CNV_status_file.txt', 'w')
        tempfile.write('gene\tCNV_status\n')
        for gene in CNV_status:
            tempfile.write(gene + '\t' + CNV_status[gene] + '\n')
        tempfile.close()
    
        # make temp summary file with targets and cnv status
        make_summary_table_target_sites(predicted_targets, 'Temp_CNV_status_file.txt', 'Temp_summary_targets.txt')
                
        # count the number of CNV genes for that study
        Num_cnv_genes = 0
        # open temp summary file for reading
        infile = open('Temp_summary_targets.txt', 'r')
        # skip header
        infile.readline()
        # loop over file
        for line in infile:
            line = line.rstrip()
            if line != '':
                line = line.split()
                if line [-1] == 'CNV':
                    Num_cnv_genes += 1
        #close file after reading
        infile.close()
        print('Num CNV genes', Num_cnv_genes)
    
        # check that study includes minimum number of cnv genes
        if Num_cnv_genes >= minimum_cnv:
            # parse the summary table into a list
            regulation = compare_miRNA_regulation('Temp_summary_targets.txt')
    
            # write regulation to file
            newfile.write(study + '\t')
            newfile.write('\t'.join(list(map(Gstr, regulation))) + '\n')
    
            print('done writing regulation for {0}'.format(study))

# close file after writing
newfile.close()
