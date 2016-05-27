# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 12:54:38 2015

@author: RJovelin
"""


# usage python3 make_summary_tables_predicted_sites.py [3UTR/5UTR/CDS] [targetscan/miranda] [True/False] [long_CNVs/all_CNVs]

# write to file tables with number of target sites, sequence length, and CNV status

from CNV_miRNAs import *
import os
import sys

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

# make a dictionary of species names : species code
species_names = {'H_sapiens': 'Hsa',  'P_troglodytes': 'Ptr', 'M_mulatta': 'Mmul',
                 'M_musculus': 'Mmus', 'R_norvegicus': 'Rno', 'B_taurus': 'Bta',
                 'C_familiaris': 'Cfa', 'G_gallus': 'Gga'}

# loop over species names
for species in species_names:
    print(species)
    
    # get the seq input file
    seq_input_file = species + '_' + domain + '_' + chromos + '_targetscan.txt'
    print(seq_input_file)
    # get the targetscan output file
    predicted_targets = species + '_' + domain + '_' + chromos + '_predicted_sites_' + predictor + '.txt'
    print(predicted_targets)
    
    # parse the predictor outputfile to get a dict {gene: [targets, seq_length, normalized_targets]}
    if predictor == 'targetscan':
        targets = parse_targetscan_output(seq_input_file, predicted_targets, 'all')
    elif predictor == 'miranda':
        targets = parse_miranda_output(seq_input_file, predicted_targets, 'all')
    print('targets', len(targets))    
    
     
    # check if species is human
    if species == 'H_sapiens':
        # make a list of files with cnv genes
        cnv_files = [species + '_GRCh37_2013-05_'  + cnv_length + '_' + chromos + '.txt',
                     species + '_GRCh37_2013-07_'  + cnv_length + '_' + chromos + '.txt',
                     species + '_GRCh37_2014_'  + cnv_length + '_' + chromos + '.txt',
                     species + '_GRCh37_2015_'  + cnv_length + '_' + chromos + '.txt']
        
        # loop over CNV gene file
        for filename in cnv_files:
            print(filename)
            # get outputfile
            if '2013-05' in filename:
                outputfile = species + '_' + domain + '_summary_' + predictor + '_' + chromos + '_' + cnv_length + '_GRCh37_2013-05.txt'
            elif '2013-07' in filename:
                outputfile = species + '_' + domain + '_summary_' + predictor + '_' + chromos + '_' + cnv_length + '_GRCh37_2013-07.txt'
            elif '2014' in filename:
                outputfile = species + '_' + domain + '_summary_' + predictor + '_' + chromos + '_' + cnv_length + '_GRCh37_2014.txt'
            elif '2015' in filename:
                outputfile = species + '_' + domain + '_summary_' + predictor + '_' + chromos + '_' + cnv_length + '_GRCh37_2015.txt'
            print(outputfile)
            
            # use this function to write summary table with target site counts and CNV status
            make_summary_table_target_sites(targets, filename, outputfile)
            
    else:
        # get CNV file
        CNV_file = species + '_' + cnv_length + '_' + chromos + '.txt' 
        print(CNV_file)
        # get outputfile
        outputfile = species + '_' + domain + '_summary_' + predictor + '_' + chromos + '_' + cnv_length + '.txt'
        print(outputfile)
    
        # use this function to write summary table with target site counts and CNV status
        make_summary_table_target_sites(targets, CNV_file, outputfile)
    
    print('done writing summary file for {0}'.format(species))
