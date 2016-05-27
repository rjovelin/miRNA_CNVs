# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 11:01:55 2015

@author: RJovelin
"""


from CNV_miRNAs import *
import os
import sys
import numpy as np
from scipy import stats
import math


# usage python3 interspecies_comparison_targets_short3UTR_genes.py options
# [5UTR/CDS]
# [targetscan/miranda]
# [True/False]
# [long_CNVs/all_CNVs]
# [shared_miRs/all_miRs]
# human_CNV_file 

# get the region to consider to predict target sites [5UTr or CDS]
domain = sys.argv[1]
print(domain)

# get the predcitor used to predictor target sites [targetscan or miranda]
predictor = sys.argv[2]
print(predictor)

# get the option to keep genes on all chromos (False) or only on assembled 
# nuclear chromosomes only from the command
keep_valid_chromos = sys.argv[3]

# check if all chromos (including unplaced, unlocated, and MT) are used
# or if only valid chromos are used 
# note: only CNVs on valid chromos are reported in DGV, so if all chromos are
# used it may introduce a bias by calling non CNV genes genes that cannot be detected

if keep_valid_chromos == 'True':
    keep_valid_chromos = True
    chromos = 'valid_chromos'
elif keep_valid_chromos == 'False':
    keep_valid_chromos = False
    chromos = 'all_chromos'
print(keep_valid_chromos, chromos)


# get the type of CNVs to consider from the command
long_CNV = sys.argv[4]
if long_CNV == 'all_CNVs':
    cnv_length = 'CNV_all_length'
elif long_CNV == 'long_CNVs':
    cnv_length = 'CNV_greater_1Kb'
print(long_CNV, cnv_length)

# get the option to consider all miRNAs or conserved miRNAs between species pairs
# [shared_miRs or all_miRs]
conservation = sys.argv[5]
print(conservation)

# get the human CNV file from the command, so that comparisons can be made 
# between all species and different version of the human DVG
human_CNV_file = sys.argv[6]
print(human_CNV_file)

# get realease version of the DGV
release_version = human_CNV_file[human_CNV_file.index('GRCh37') : human_CNV_file.index('_CNV')]

# get outputfile
outputfile = 'Interspecies_comp_targets_' + predictor + '_' + domain + '_short3UTR_' + chromos + '_' + cnv_length + '_' + conservation + '_' + release_version + '.txt'
print(outputfile)

# open file for writing
newfile = open(outputfile, 'w')

# write header to file
newfile.write('Species1' + '\t' + 'Species2' + '\t')

# write number of miRNAs if using shared miRNAs
if conservation  == 'shared_miRs':
    newfile.write('shared_miRNAs' + '\t')

newfile.write('\t'.join(['N_CNV_genes_sp1', 'N_CNV_genes_sp2',
                         'CNV_mean_targets_sp1', 'CNV_SEM_targets_sp1',
                         'CNV_mean_targets_sp2', 'CNV_SEM_targets_sp2', 'P_diff_CNV_targets',
                         'N_nonCNV_genes_sp1', 'N_nonCNV_genes_sp2', 
                         'nonCNV_mean_targets_sp1', 'nonCNV_SEM_targets_sp1',
                         'nonCNV_mean_targets_sp2', 'nonCNV_SEM_targets_sp2', 'P_diff_nonCNV_targets',
                         'CNV_mean_normalized_targets_sp1', 'CNV_SEM_normalized_targets_sp1',
                         'CNV_mean_normalized_targets_sp2', 'CNV_SEM_normalized_targets_sp2', 'P_diff_norm_CNV_targets',
                         'nonCNV_mean_normalized_targets_sp1', 'nonCNV_SEM_normalized_targets_sp1',
                         'nonCNV_mean_normalized_targets_sp2', 'nonCNV_SEM_normalized_targets_sp2', 'P_diff_norm_nonCNV_targets',
                         'CNV_mean_seq_length_sp1', 'CNV_SEM_seq_length_sp1',
                         'CNV_mean_seq_length_sp2', 'CNV_SEM_seq_length_sp2', 'P_diff_seq_CNV_length',
                         'nonCNV_mean_seq_length_sp1', 'nonCNV_SEM_seq_length_sp1',
                         'nonCNV_mean_seq_length_sp2', 'nonCNV_SEM_seq_length_sp2', 'P_diff_seq_nonCNV_length']) + '\n')                        

# create a lamda function to transform value into string
Gstr = lambda x: str(x)

# make a list of species names
species = ['H_sapiens', 'P_troglodytes', 'M_mulatta', 'M_musculus', 'R_norvegicus', 'B_taurus', 'C_familiaris', 'G_gallus']


# loop over species
for i in range(len(species) -1):
    # compare pairs of species
    for j in range(i+1, len(species)):
        print(species[i], species[j])
        # get genome files        
        genome_sp1 = species[i] + '_genome.txt'
        genome_sp2 = species[j] + '_genome.txt'
        # get chromos files
        valid_chromos_sp1 = species[i] + '_valid_chromos.txt'
        valid_chromos_sp2 = species[j] + '_valid_chromos.txt'
        # get GFF annotation files
        GFF_sp1 = species[i] + '.gff3'
        GFF_sp2 = species[j] + '.gff3'
        print(genome_sp1, genome_sp2)
        print(valid_chromos_sp1, valid_chromos_sp2)
        print(GFF_sp1, GFF_sp2)
        
        # get targetscan sequence input file
        seq_input_sp1 = species[i] + '_' + domain + '_' + chromos + '_targetscan.txt'
        seq_input_sp2 = species[j] + '_' + domain + '_' + chromos + '_targetscan.txt'
        print(seq_input_sp1, seq_input_sp2)
        
        # get predictor output
        predicted_targets_sp1 = species[i] + '_' + domain + '_' + chromos + '_predicted_sites_' + predictor + '.txt'
        predicted_targets_sp2 = species[j] + '_' + domain + '_' + chromos + '_predicted_sites_' + predictor + '.txt'
        print(predicted_targets_sp1, predicted_targets_sp2)        
                
        # get UTR files
        UTR_sp1 = species[i] + '_3UTR_length_' + chromos + '.txt'      
        UTR_sp2 = species[j] + '_3UTR_length_' + chromos + '.txt'
        print(UTR_sp1, UTR_sp2)
         
        # get CNV files
        if species[i] == 'H_sapiens':
            CNV_file_sp1 = human_CNV_file
            CNV_file_sp2 = species[j] + '_' + cnv_length + '_' + chromos + '.txt'
        else:
            CNV_file_sp1 = species[i] + '_' + cnv_length + '_' + chromos + '.txt'
            CNV_file_sp2 = species[j] + '_' + cnv_length + '_' + chromos + '.txt'
        print(CNV_file_sp1, CNV_file_sp2)
        
        # get mature mirna files
        mature_sp1 = species[i] + '_mature.txt'
        mature_sp2 = species[j] + '_mature.txt'        
        print(mature_sp1, mature_sp2)
        
        # get a set of shared seeds between species pair
        shared_seeds = find_conserved_mirna_families(mature_sp1, mature_sp2)
        
        # get CNV gene status
        CNV_status1 = sort_genes_CNV_status(CNV_file_sp1)
        CNV_status2 = sort_genes_CNV_status(CNV_file_sp2)
        print('CNV status', len(CNV_status1), len(CNV_status2))
        
        # parse predictor outputfile {gene: [N_sites, seq_length, N_sites/seq_length]}
        # check predictor
        if predictor == 'targetscan':
            # check if all mirna families are considered or only families shared by the 2 species
            if conservation == 'shared_miRs':
                # filter targets for shared mirnas
                targets_sp1 = parse_targetscan_output(seq_input_sp1, predicted_targets_sp1, 'conserved', shared_seeds, mature_sp1)
                targets_sp2 = parse_targetscan_output(seq_input_sp2, predicted_targets_sp2, 'conserved', shared_seeds, mature_sp2)
            elif conservation == 'all_miRs':
                # take all sites for all mirnas
                targets_sp1 = parse_targetscan_output(seq_input_sp1, predicted_targets_sp1, 'all')
                targets_sp2 = parse_targetscan_output(seq_input_sp2, predicted_targets_sp2, 'all')
        
        elif predictor == 'miranda':
            # check if all mirna families are considered or only families shared by the 2 species
            if conservation == 'shared_miRs':
                targets_sp1 = parse_miranda_output(seq_input_sp1, predicted_targets_sp1, 'conserved', shared_seeds, mature_sp1)
                targets_sp2 = parse_miranda_output(seq_input_sp2, predicted_targets_sp2, 'conserved', shared_seeds, mature_sp2)
            elif conservation == 'all_miRs':
                targets_sp1 = parse_miranda_output(seq_input_sp1, predicted_targets_sp1, 'all')
                targets_sp2 = parse_miranda_output(seq_input_sp2, predicted_targets_sp2, 'all')
        print('targets', len(targets_sp1), len(targets_sp2))
        
        # sort genes by UTR length
        UTR_length_sp1 = sort_genes_3UTR_length(UTR_sp1)
        UTR_length_sp2 = sort_genes_3UTR_length(UTR_sp2)
        print('UTR length', len(UTR_length_sp1), len(UTR_length_sp2))
        
        # create list of number of targets, normalized targets, sequence length for genes in CNV
        Sp1_CNV_sites, Sp2_CNV_sites, Sp1_CNV_norm_sites, Sp2_CNV_norm_sites, Sp1_CNV_length, Sp2_CNV_length = [], [], [], [], [], []
        
        # create list of number of targets, normalized targets, sequence length for genes in non CNV
        Sp1_nonCNV_sites, Sp2_nonCNV_sites, Sp1_nonCNV_norm_sites, Sp2_nonCNV_norm_sites, Sp1_nonCNV_length, Sp2_nonCNV_length = [], [], [], [], [], []        
        
        # loop over genes in targets sp1
        for gene in targets_sp1:
            # record only genes with short 3'UTR
            if UTR_length_sp1[gene] == 'short':
                # check if gene in CNV or not
                if CNV_status1[gene] == 'CNV':
                    Sp1_CNV_sites.append(targets_sp1[gene][0])
                    Sp1_CNV_norm_sites.append(targets_sp1[gene][2])
                    Sp1_CNV_length.append(targets_sp1[gene][1])
                elif CNV_status1[gene] == 'not_CNV':
                    Sp1_nonCNV_sites.append(targets_sp1[gene][0])
                    Sp1_nonCNV_norm_sites.append(targets_sp1[gene][2])
                    Sp1_nonCNV_length.append(targets_sp1[gene][1])
        
        # loop over genes in targets sp2
        for gene in targets_sp2:
            # record only genes with short 3'UTR
            if UTR_length_sp2[gene] == 'short':
                # check if gene in CNV or not
                if CNV_status2[gene] == 'CNV':
                    Sp2_CNV_sites.append(targets_sp2[gene][0])
                    Sp2_CNV_norm_sites.append(targets_sp2[gene][2])
                    Sp2_CNV_length.append(targets_sp2[gene][1])
                elif CNV_status2[gene] == 'not_CNV':
                    Sp2_nonCNV_sites.append(targets_sp2[gene][0])
                    Sp2_nonCNV_norm_sites.append(targets_sp2[gene][2])
                    Sp2_nonCNV_length.append(targets_sp2[gene][1])
                
                    
        print('CNV targets', len(Sp1_CNV_sites), len(Sp2_CNV_sites))
        print('non CNV targets', len(Sp1_nonCNV_sites), len(Sp2_nonCNV_sites))
        
        
        # write content to file
        newfile.write(species[i] + '\t' + species[j] + '\t')
        
        # write number of miRNAs if using shared miRNAs
        if conservation == 'shared_miRs':
            newfile.write(str(len(shared_seeds)) + '\t')
        
        newfile.write(str(len(Sp1_CNV_sites)) + '\t' + str(len(Sp2_CNV_sites)) + '\t')
        newfile.write(str(np.mean(Sp1_CNV_sites)) + '\t' + str(np.std(Sp1_CNV_sites) / math.sqrt(len(Sp1_CNV_sites))) + '\t')
        newfile.write(str(np.mean(Sp2_CNV_sites)) + '\t' + str(np.std(Sp2_CNV_sites) / math.sqrt(len(Sp2_CNV_sites))) + '\t')
        newfile.write(str(stats.ranksums(Sp1_CNV_sites, Sp2_CNV_sites)[1]) + '\t')
        newfile.write(str(len(Sp1_nonCNV_sites)) + '\t' + str(len(Sp2_nonCNV_sites)) + '\t')
        newfile.write(str(np.mean(Sp1_nonCNV_sites)) + '\t' + str(np.std(Sp1_nonCNV_sites) / math.sqrt(len(Sp1_nonCNV_sites))) + '\t')
        newfile.write(str(np.mean(Sp2_nonCNV_sites)) + '\t' + str(np.std(Sp2_nonCNV_sites) / math.sqrt(len(Sp2_nonCNV_sites))) + '\t')
        newfile.write(str(stats.ranksums(Sp1_nonCNV_sites, Sp2_nonCNV_sites)[1]) + '\t')
        newfile.write(str(np.mean(Sp1_CNV_norm_sites)) + '\t' + str(np.std(Sp1_CNV_norm_sites) / math.sqrt(len(Sp1_CNV_norm_sites))) + '\t')
        newfile.write(str(np.mean(Sp2_CNV_norm_sites)) + '\t' + str(np.std(Sp2_CNV_norm_sites) / math.sqrt(len(Sp2_CNV_norm_sites))) + '\t')
        newfile.write(str(stats.ranksums(Sp1_CNV_norm_sites, Sp2_CNV_norm_sites)[1]) + '\t')
        newfile.write(str(np.mean(Sp1_nonCNV_norm_sites)) + '\t' + str(np.std(Sp1_nonCNV_norm_sites) / math.sqrt(len(Sp1_nonCNV_norm_sites))) + '\t')
        newfile.write(str(np.mean(Sp2_nonCNV_norm_sites)) + '\t' + str(np.std(Sp2_nonCNV_norm_sites) / math.sqrt(len(Sp2_nonCNV_norm_sites))) + '\t')
        newfile.write(str(stats.ranksums(Sp1_nonCNV_norm_sites, Sp2_nonCNV_norm_sites)[1]) + '\t')              
        newfile.write(str(np.mean(Sp1_CNV_length)) + '\t' + str(np.std(Sp1_CNV_length) / math.sqrt(len(Sp1_CNV_length))) + '\t')
        newfile.write(str(np.mean(Sp2_CNV_length)) + '\t' + str(np.std(Sp2_CNV_length) / math.sqrt(len(Sp2_CNV_length))) + '\t')
        newfile.write(str(stats.ranksums(Sp1_CNV_length, Sp2_CNV_length)[1]) + '\t')
        newfile.write(str(np.mean(Sp1_nonCNV_length)) + '\t' + str(np.std(Sp1_nonCNV_length) / math.sqrt(len(Sp1_nonCNV_length))) + '\t')
        newfile.write(str(np.mean(Sp2_nonCNV_length)) + '\t' + str(np.std(Sp2_nonCNV_length) / math.sqrt(len(Sp2_nonCNV_length))) + '\t')
        newfile.write(str(stats.ranksums(Sp1_nonCNV_length, Sp2_nonCNV_length)[1]) + '\n')
        
        print('done writing regulation for {0} and {1}'.format(species[i], species[j]))

# close file after writing
newfile.close()