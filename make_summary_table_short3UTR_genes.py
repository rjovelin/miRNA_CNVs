# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 22:59:58 2015

@author: Richard
"""

# make summary tables of target sites predicted in CDS or 5' UTR for
# genes with short 3' UTRs using targetscan or miranda as predictor


# usage python3 make_summary_table_short3UTR_genes.py
# [5UTR/CDS]
# [targetscan/miranda]
# [True/False]
# [long_CNVs/all_CNVs]


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
        # use DGV 2015 release 
        CNV_file = 'H_sapiens_GRCh37_2015_CNV_all_length_valid_chromos.txt'
        outputfile = species + '_short3UTRgenes_' + domain + '_summary_' + predictor + '_' + chromos + '_' + cnv_length + '_GRCh37_2015.txt'
    else:
        # get CNV file
        CNV_file = species + '_' + cnv_length + '_' + chromos + '.txt' 
        # get outputfile
        outputfile = species + '_short3UTRgenes_' + domain + '_summary_' + predictor + '_' + chromos + '_' + cnv_length + '.txt'
    
    print(CNV_file)
    # get file with UTR status (short, long)
    UTR_file = species + '_3UTR_length_' + chromos + '.txt'
    print(UTR_file)
    # sort genes by UTR length
    UTR_length = sort_genes_3UTR_length(UTR_file)
    print('UTR length', len(UTR_length))
    
    # get CNV gene status
    CNV_status = sort_genes_CNV_status(CNV_file)
    print('CNV status', len(CNV_status))
    
    # write number of target sites, sequence length, normalized targets and CNV status to file
    newfile = open(outputfile, 'w')
    # write header
    newfile.write('\t'.join(['Gene', 'N_targets', 'Sequence_length', 'N_targets_normalized', 'CNV_status']) + '\n')
    
    # loop over genes in targets
    for gene in targets:
        # record only genes with short 3'UTR
        if gene in UTR_length and UTR_length[gene] == 'short':
            # write number of sites, sequence length and normalized number of targets to file
            newfile.write('\t'.join([gene, str(targets[gene][0]), str(targets[gene][1]), str(targets[gene][2])]) + '\t')
            # write CNV status
            newfile.write(CNV_status[gene] + '\n')
    # close file after writing
    newfile.close()
    print('done writing targets to file')
    
