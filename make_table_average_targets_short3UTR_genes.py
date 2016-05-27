# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 14:17:55 2015

@author: RJovelin
"""

from CNV_miRNAs import *
import os
import sys
import numpy as np
from scipy import stats
import math


# usage python3 make_table_average_targets_short3UTR_genes.py [5UTR/CDS] [targetscan/miranda] [True/False] [long_CNVs/all_CNVs]

# get the region to consider to predict target sites [5UTr or CDS]
domain = sys.argv[1]
print(domain)

# get predictor from command [targetscan or miranda]
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
outputfile = 'Average_targets_' + domain + '_short3UTR_' + chromos + '_' + cnv_length + '_' + predictor + '.txt'
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
    
    # get genome file
    genome_file = species + '_genome.txt'
    # get GFF annotation file
    GFF_file = species + '.gff3'
    # get valid chromos file
    valid_chromos = species + '_valid_chromos.txt'
    print(genome_file)
    print(GFF_file)
    print(valid_chromos)
    
    # get targetscan sequence input file
    seq_input_file = species + '_' + domain + '_' + chromos + '_targetscan.txt'
    print(seq_input_file)    
    # get outputfile with predicted sites
    predicted_target_file = species + '_' + domain + '_' + chromos + '_predicted_sites_' + predictor + '.txt'
    print(predicted_target_file)
    # get file with UTR status (short, long)
    UTR_file = species + '_3UTR_length_' + chromos + '.txt'
    print(UTR_file)
    
    # create a dict to store the number of targets {gene: [N_sites, seq_length, N_sites/seq_length]}
    if predictor == 'targetscan':
        targets = parse_targetscan_output(seq_input_file, predicted_target_file, 'all')
    elif predictor == 'miranda':
        targets = parse_miranda_output(seq_input_file, predicted_target_file, 'all')
    print('targets', len(targets))
    
    
    # sort genes by UTR length
    UTR_length = sort_genes_3UTR_length(UTR_file)
    print('UTR length', len(UTR_length))
    
    # check if species is human
    if species == 'H_sapiens':
        # make a list of files with cnv genes
        cnv_files = [species + '_GRCh37_2013-05_'  + cnv_length + '_' + chromos + '.txt',
                     species + '_GRCh37_2013-07_'  + cnv_length + '_' + chromos + '.txt',
                     species + '_GRCh37_2014_'  + cnv_length + '_' + chromos + '.txt',
                     species + '_GRCh37_2015_'  + cnv_length + '_' + chromos + '.txt']
        
        # compute average target sites for short 3'UTR genes using each version of the DVG
        for filename in cnv_files:
            print(filename)
            
            # get CNV gene status
            CNV_status = sort_genes_CNV_status(filename)
            print('CNV status', len(CNV_status))
            
            # get DGV version
            DGV_version = filename[filename.index('GRCh'): filename.index('_CNV_')]            
            print(DGV_version)
            
            # write a temporary file with number of target sites, sequence length, and CNV status
            tempfile = open('Temp_file_targets.txt', 'w')
            # write header
            tempfile.write('\t'.join(['Gene', 'N_targets', 'Sequence_length', 'N_targets_normalized', 'CNV_status']) + '\n')
            
            # loop over genes in targets
            for gene in targets:
                # record only genes with short 3'UTR
                if gene in UTR_length and UTR_length[gene] == 'short':
                    # write number of sites, sequence length and normalized number of targets to file
                    tempfile.write('\t'.join([gene, str(targets[gene][0]), str(targets[gene][1]), str(targets[gene][2])]) + '\t')
                    # write CNV status
                    tempfile.write(CNV_status[gene] + '\n')
            # close file after writing
            tempfile.close()
            print('done writing targets to file')
            
            # parse the summary table into a list
            regulation = compare_miRNA_regulation('Temp_file_targets.txt')
            
            # write regulation to file
            newfile.write(species + '_' + DGV_version + '\t')
            newfile.write('\t'.join(list(map(Gstr, regulation))) + '\n')
    
            print('done writing regulation for {0}'.format(species + '_' + DGV_version))
            
            
    else:
        # get CNV file
        CNV_file = species + '_' + cnv_length + '_' + chromos + '.txt' 
        print(CNV_file)
        
        # get CNV gene status
        CNV_status = sort_genes_CNV_status(CNV_file)
        print('CNV status', len(CNV_status))
    
        # write a temporary file with number of target sites, sequence length, and CNV status
        tempfile = open('Temp_file_targets.txt', 'w')
        # write header
        tempfile.write('\t'.join(['Gene', 'N_targets', 'Sequence_length', 'N_targets_normalized', 'CNV_status']) + '\n')
    
        # loop over genes in targets
        for gene in targets:
            # record only genes with short 3'UTR
            if gene in UTR_length and UTR_length[gene] == 'short':
                # write number of sites, sequence length and normalized number of targets to file
                tempfile.write('\t'.join([gene, str(targets[gene][0]), str(targets[gene][1]), str(targets[gene][2])]) + '\t')
                # write CNV status
                tempfile.write(CNV_status[gene] + '\n')
        # close file after writing
        tempfile.close()
        print('done writing targets to file')

        # parse the summary table into a list
        regulation = compare_miRNA_regulation('Temp_file_targets.txt')
    
        # write regulation to file
        newfile.write(species + '\t')
        newfile.write('\t'.join(list(map(Gstr, regulation))) + '\n')
    
        print('done writing regulation for {0}'.format(species))

# close file after writing
newfile.close()


