# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 14:38:27 2016

@author: RJovelin
"""

# use this script to compute the correlation between the number of mirna target sites
# in regions of genes with short 3'UTRs and the number of miRNAs per species

# use Agg backend on server without X server
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rc
rc('mathtext', default='regular')
# import modules
import numpy as np
from scipy import stats
import math
import os
import sys
import copy
# import custom modules
from CNV_miRNAs import *


# usage python3 CorrelationTargetSitesShort3UTRs.py options
# ['spearman', 'pearson']: use spearman or pearson correlation
# ['mature', 'family', 'unique']: perform correlation on the number of mature MIRNAs, miRNA families or unique families in each species
# [7/15]: 3'UTR length, genes with 3'UTR length < 7bp or < 15bp are considered short 3'UTR genes

# get the correlation type from command
correlation = sys.argv[1]
assert correlation in ['spearman', 'pearson'], 'correlation should be spearman or pearson'
# get the variable to perform correlation
miRvar = sys.argv[2]
assert miRvar in ['mature', 'family', 'unique'], 'miRvar should be mature, family or unique'
# get the minimum 3'UTR length
L = int(sys.argv[3])
assert L in [15, 7], 'minimum 3UTR length is not correct'

# use all chromos (including unplaced, unlocated, and MT) or only valid chromos 
# note: only CNVs on valid chromos are reported in DGV, so if all chromos are
# used it may introduce a bias by calling non CNV genes genes that cannot be detected

# keep genes on assembled nuclear chromosomes
chromos = 'valid_chromos'
# consider all CNVs
cnv_length = 'CNV_all_length'

# get the number of mirna per species {species: [N_mature, N_families, N_unique_families]}
miRNASpecies = {}
infile = open('mirna_family_species.txt')
infile.readline()
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split()
        miRNASpecies[line[0].capitalize()] = list(map(lambda x: int(x), line[1:]))
infile.close()
print('recorded the number of miRNAs in each speacies')


# make a dictionary of species names : species code
species_codes = {'H_sapiens': 'Hsa', 'P_troglodytes': 'Ptr', 'M_mulatta': 'Mml',
                 'M_musculus': 'Mmu', 'B_taurus': 'Bta', 'G_gallus':'Gga'}

# make a dictionnary {species: {domain: mean number of normalized targets}}
# for each predictor, targetscan and miranda
targetscan, miranda = {}, {}

# make lists of predictors and regions
predictors = ['targetscan', 'miranda']
regions = ['5UTR', 'CDS']
# make a list of species names to loop from
SpeciesNames = ['H_sapiens', 'P_troglodytes', 'M_mulatta', 'M_musculus', 'B_taurus', 'G_gallus']


# loop over species names
for species in species_codes:
    # loop over predictors
    for i in range(len(predictors)):
        # loop over domain
        for domain in regions:
            # get the seq input file
            seq_input_file = species + '_' + domain + '_' + chromos + '_targetscan.txt'
            # get the predicted targets output file
            predicted_targets = species + '_' + domain + '_' + chromos + '_predicted_sites_' + predictors[i] + '.txt'
            # parse the predictor outputfile to get a dict {gene: [targets, seq_length, normalized_targets]}
            if i == 0:
                targets = parse_targetscan_output(seq_input_file, predicted_targets, 'all')
            elif i == 1:
                targets = parse_miranda_output(seq_input_file, predicted_targets, 'all')
            # get file with UTR status (short, long)
            UTR_file = species + '_3UTR_length_' + chromos + '.txt'
            # create a dict {gene : "short" (or "long")}
            UTR_length = sort_genes_3UTR_length(UTR_file, L)
            # make a list of targets per nucleotide per gene
            TargetNum = []
            for gene in targets:
                # record gene with short 3'UTR
                if gene in UTR_length and UTR_length[gene] == 'short':
                    TargetNum.append(targets[gene][2])
            # initialize inner dict
            if i == 0:
                # check if species in dict
                if species_codes[species] not in targetscan:
                    targetscan[species_codes[species]] = {}
                # populate dict with domain mean targets pairs
                targetscan[species_codes[species]][domain] = np.mean(TargetNum)
            elif i == 1:
                # check if species in dict
                if species_codes[species] not in miranda:
                    miranda[species_codes[species]] = {}
                # populate dict with domain mean target pairs
                miranda[species_codes[species]][domain] = np.mean(TargetNum)
print('recorded the mean number of targets')                

# create parallel lists with mean number of targets and number of mirnas in each species
TargetscanUTR, TargetscanCDS, MirandaUTR, MirandaCDS, mirnas = [], [], [], [], []
for species in miRNASpecies:
    if miRvar == 'mature':
        mirnas.append(miRNASpecies[species][0])
    elif miRvar == 'family':
        mirnas.append(miRNASpecies[species][1])
    elif miRvar == 'unique':
        mirnas.append(miRNASpecies[species][2])
    TargetscanUTR.append(targetscan[species]['5UTR'])
    MirandaUTR.append(mirnada[species]['5UTR'])
    TargetscanCDS.append(targetscan[species]['CDS'])
    MirandaCDS.append(miranda[species]['CDS'])
        
# perform correlations between mean number of targets and mirna numbers
if correlation == 'spearman':
    CorrelTargetscanUTR = stats.spearmanr(TargetscanUTR, mirnas)[0]
    CorrelMirandaUTR = stats.spearmanr(MirandaUTR, mirnas)[0]
    CorrelTargetscanCDS = stats.spearmanr(TargetscanCDS, mirnas)[0]
    CorrelMirandaCDS = stats.spearmanr(MirandaCDS, mirnas)[0]
elif correlation == 'pearson':
    CorrelTargetscanUTR = stats.pearsonr(TargetscanUTR, mirnas)[0]
    CorrelMirandaUTR = stats.pearsonr(MirandaUTR, mirnas)[0]
    CorrelTargetscanCDS = stats.pearsonr(TargetscanCDS, mirnas)[0]
    CorrelMirandaCDS = stats.pearsonr(MirandaCDS, mirnas)[0]
print('performed correlations') 

 
# print results to screen
if miRvar == 'mature':
    print('correlations between mean number of targets and number of mature miRNAs') 
elif miRvar == 'family':
    print('correlations between mean number of targets and number of miRNA families')
elif miRvar == 'unique':
    print('correlations between mean number of targets and number of unique miRNA families')
print('\n')

print('5UTR', 'targetscan', CorrelTargetscanUTR)
print('5UTR', 'miranda', CorrelMirandaUTR)
print('CDS', 'targetscan', CorrelTargetscanCDS)
print('CDS', 'miranda', CorrelMirandaCDS)