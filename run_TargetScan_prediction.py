# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 15:31:41 2015

@author: RJovelin
"""

import os


# usage python3 run_TargetScan_prediction.py [3UTR/5UTR/CDS] [True/False]


import os
import sys

# get the region to consider to predict target sites [UTR or 5UTr or CDS]
domain = sys.argv[1]


# get the option to keep genes on all chromos (False) or only on assembled 
# nuclear chromosomes only from the command
keep_valid_chromos = sys.argv[2]

# check if all chromos (including unplaced, unlocated, and MT) are used
# or if only valid chromos are used 
# note: only CNVs on valid chromos are reported in DGV, so if all chromos are
# used it may introduce a bias by calling non CNV genes genes that cannot be detected

if keep_valid_chromos == 'True':
    # get the file and the end of the targetscan input file name
    input_file_end = '_valid_chromos_targetscan.txt'
elif keep_valid_chromos == 'False':
    input_file_end = '_all_chromos_targetscan.txt'
    
# make a dictionary of species names : species code
species_names = {'H_sapiens': 'Hsa',  'P_troglodytes': 'Ptr', 'M_mulatta': 'Mmul',
                 'M_musculus': 'Mmus', 'R_norvegicus': 'Rno', 'B_taurus': 'Bta',
                 'C_familiaris': 'Cfa', 'G_gallus': 'Gga'}

# loop over species
for species in species_names:
    print(species)
    
    # get the miRNA input file
    miRNA_input = species + '_miRFam_targetscan.txt'
    print(miRNA_input)
    
    # get the sequence input file
    sequence_input_file = species + '_' + domain + input_file_end
    print(sequence_input_file)
    
    # get outputfile
    # check if which chromosomes were considered
    if keep_valid_chromos == 'True':
        outputfile = species + '_' + domain + '_valid_chromos_predicted_sites_targetscan.txt'
    elif keep_valid_chromos == 'False':
        outputfile = species + '_' + domain + '_all_chromos_predicted_sites_targetscan.txt'
    print(outputfile)
    
    # build TargetScan command 
    os.system('perl targetscan_60.pl ' + miRNA_input + ' ' + sequence_input_file + ' ' + outputfile)
    
    