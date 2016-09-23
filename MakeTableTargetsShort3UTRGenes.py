# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 15:27:19 2016

@author: RJovelin
"""

# use this script to generate a table with target sites in CNV and non-CNv genes 
# that have short 3'UTRs 

# usage MakeTableTargetShort3UTRGenes.py [options]
# [7/15]: 3'UTR length, genes with 3'UTR length < 7bp or < 15bp are considered short 3'UTR genes

# import modules
import numpy as np
from scipy import stats
import math
import os
import sys
import copy
# import custom modules
from CNV_miRNAs import *


# get the minimum 3'UTR length
L = int(sys.argv[1])
assert L in [15, 7], 'minimum 3UTR length is not correct'
print(L)

# use all chromos (including unplaced, unlocated, and MT) or only valid chromos 
# note: only CNVs on valid chromos are reported in DGV, so if all chromos are
# used it may introduce a bias by calling non CNV genes genes that cannot be detected

# keep genes on assembled nuclear chromosomes
chromos = 'valid_chromos'
# consider all CNVs
cnv_length = 'CNV_all_length'

# make a dictionary of species names : species code
species_codes = {'H_sapiens': 'Hsa', 'P_troglodytes': 'Ptr', 'M_mulatta': 'Mml',
                 'M_musculus': 'Mmu', 'B_taurus': 'Bta', 'G_gallus':'Gga'}

# make a dictionnary {species: {domain: {gene: [N_targets, Sequence_length, N_targets_normalized, CNV_status]}}}
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
            print(i, predictors[i], species, domain)
            # get the seq input file
            seq_input_file = species + '_' + domain + '_' + chromos + '_targetscan.txt'
            print(seq_input_file)
            # get the predicted targets output file
            predicted_targets = species + '_' + domain + '_' + chromos + '_predicted_sites_' + predictors[i] + '.txt'
            print(predicted_targets)
            # parse the predictor outputfile to get a dict {gene: [targets, seq_length, normalized_targets]}
            if i == 0:
                targets = parse_targetscan_output(seq_input_file, predicted_targets, 'all')
            elif i == 1:
                targets = parse_miranda_output(seq_input_file, predicted_targets, 'all')
            print('targets', len(targets))    
            # check if species is human
            if species == 'H_sapiens':
                # use DGV 2015 release 
                CNV_file = 'H_sapiens_GRCh37_2015_CNV_all_length_valid_chromos.txt'
            else:
                CNV_file = species + '_' + cnv_length + '_' + chromos + '.txt' 
            print(CNV_file)
            # get CNV gene status
            CNV_status = sort_genes_CNV_status(CNV_file)
            print('CNV status', len(CNV_status))
            # add CNV status
            for gene in targets:
                targets[gene].append(CNV_status[gene])
            # get file with UTR status (short, long)
            UTR_file = species + '_3UTR_length_' + chromos + '.txt'
            print(UTR_file)
            # create a dict {gene : "short" (or "long")}
            UTR_length = sort_genes_3UTR_length(UTR_file, L)
            print('UTR length', len(UTR_length))
            # initialize inner dict
            if i == 0:
                # check if species in dict
                if species not in targetscan:
                    targetscan[species] = {}
                #initialize inner dict with domain as key
                targetscan[species][domain] = {}
                for gene in targets:
                    # record gene with short 3'UTR
                    if gene in UTR_length and UTR_length[gene] == 'short':
                         targetscan[species][domain][gene] = copy.deepcopy(targets[gene])
            elif i == 1:
                # check if species in dict
                if species not in miranda:
                    miranda[species] = {}
                # initialize inner dict with domain as key
                miranda[species][domain] = {}
                for gene in targets:
                    # record gene with short 3'UTR
                    if gene in UTR_length and UTR_length[gene] == 'short':
                        miranda[species][domain][gene] = copy.deepcopy(targets[gene])

# check that the same number of species are recorded for miranda and targetscan
assert len(targetscan) == len(miranda), 'different number of species depending on predictor'
# check that list values have the correct number of items
for species in targetscan:
    for domain in targetscan[species]:
        for gene in targetscan[species][domain]:
            assert len(targetscan[species][domain][gene]) == 4, 'gene in targetscan does not have all required values'
        for gene  in miranda[species][domain]:
            assert len(miranda[species][domain][gene]) == 4, 'gene in miranda does not have all required values'


# create a dict of {species name : {domain: [CNV], [non-CNV]]}}
DataTargetscan, DataMiranda = {}, {}

# loop over species names, get the number of normalized sites for CNV and non-CNV genes
for species in targetscan:
    # initialise inner dict
    DataTargetscan[species] = {}
    # loop over domain
    for domain in targetscan[species]:
        # initialize list value
        DataTargetscan[species][domain] = [[], []]    
        # populate inner lists with number of miRNA target sites per nucleotide
        for gene in targetscan[species][domain]:
            if targetscan[species][domain][gene][-1] == 'CNV':
                DataTargetscan[species][domain][0].append(targetscan[species][domain][gene][2])
            elif targetscan[species][domain][gene][-1] == 'not_CNV':
                DataTargetscan[species][domain][1].append(targetscan[species][domain][gene][2])
for species in miranda:
    # initialize inner dict
    DataMiranda[species] = {}
    # loop over domain
    for domain in miranda[species]:
        # initialize list value
        DataMiranda[species][domain] = [[], []]
        # populate inner lists with number of mirna target sites per nucleotide
        for gene in miranda[species][domain]:
            if miranda[species][domain][gene][-1] == 'CNV':
                DataMiranda[species][domain][0].append(miranda[species][domain][gene][2])
            elif miranda[species][domain][gene][-1] == 'not_CNV':
                DataMiranda[species][domain][1].append(miranda[species][domain][gene][2])
print('generated lists of target sites for CNV and non-CNV genes')


# create a dict with gene numbers {species: {domain: [N CNV genes targetscan, N CNV genes miranda]}}
GeneNumbersTargetscan, GeneNumbersMiranda = {}, {}
for species in targetscan:
    GeneNumbersTargetscan[species] = {}
    for domain in targetscan[species]:
        GeneNumbersTargetscan[species][domain] = [0, 0]
        for gene in targetscan[species][domain]:
            if targetscan[species][domain][gene][-1] == 'CNV':
                GeneNumbersTargetscan[species][domain][0] += 1
            elif targetscan[species][domain][gene][-1] == 'not_CNV':
                GeneNumbersTargetscan[species][domain][1] += 1
for species in miranda:
    GeneNumbersMiranda[species] = {}
    for domain in miranda[species]:
        GeneNumbersMiranda[species][domain] = [0, 0]
        for gene in miranda[species][domain]:
            if miranda[species][domain][gene][-1] == 'CNV':
                GeneNumbersMiranda[species][domain][0] += 1
            elif miranda[species][domain][gene][-1] == 'not_CNV':
                GeneNumbersMiranda[species][domain][1] += 1
print('counted CNV and non-CNV genes')
            
            
# perform stattistical tests between CNV and non-CNV genes
# create a dict to store results {species: {domain: significance level}}
PvalTargetscan, PvalMiranda = {}, {}
for species in DataTargetscan:
    PvalTargetscan[species], PvalMiranda[species] = {}, {}
    for domain in DataTargetscan[species]:
        P = stats.ranksums(DataTargetscan[species][domain][0], DataTargetscan[species][domain][1])[1]
        if P >= 0.05:
            PvalTargetscan[species][domain] = ''
        elif P < 0.05 and P >= 0.01:
            PvalTargetscan[species][domain] = '*'
        elif P < 0.01 and P >= 0.001:
            PvalTargetscan[species][domain] = '**'
        elif P < 0.001:
            PvalTargetscan[species][domain] = '***'
    for domain in DataMiranda[species]:
        P = stats.ranksums(DataMiranda[species][domain][0], DataMiranda[species][domain][1])[1]
        if P >= 0.05:
            PvalMiranda[species][domain] = ''
        elif P < 0.05 and P >= 0.01:
            PvalMiranda[species][domain] = '*'
        elif P < 0.01 and P >= 0.001:
            PvalMiranda[species][domain] = '**'
        elif P < 0.001:
            PvalMiranda[species][domain] = '***'
print('compared CNV and non-CNV genes')


# write output to file
# get outputfile name
outputfile = 'TableTargetsShort3UTRgenes_' + str(L) + 'bp_' + chromos + '_' + cnv_length + '.txt'
# open file for writing
newfile = open(outputfile, 'w')
# write table header
newfile.write('\t'.join(['', '', '', '', 'TargetScan', '', '', 'miRanda', '', '']) + '\n')
newfile.write('\t'.join(['', '', 'CNV', 'Non-CNV', 'CNV', 'Non-CNV', '', 'CNV', 'Non-CNV', '']) + '\n')
newfile.write('\t'.join(['Sp', 'Domain', 'Na', 'Na', 'Mean', 'Mean', 'D (%)b', 'Mean', 'Mean', 'D (%)b']) + '\n')

# loop over domain
for domain in regions:
    # loop over species
    for species in SpeciesNames:
        # create the line to write
        line = [species_codes[species], domain, str(GeneNumbersTargetscan[species][domain][0]), str(GeneNumbersTargetscan[species][domain][1]),
                str(round(np.mean(DataTargetscan[species][domain][0]), 4)), str(round(np.mean(DataTargetscan[species][domain][1]), 4)),
                str(round((1 - np.mean(DataTargetscan[species][domain][1])/np.mean(DataTargetscan[species][domain][0])) * 100, 2)) + PvalTargetscan[species][domain],
                str(round(np.mean(DataMiranda[species][domain][0]), 4)), str(round(np.mean(DataMiranda[species][domain][1]), 4)),
                str(round((1 - np.mean(DataMiranda[species][domain][1])/np.mean(DataMiranda[species][domain][0])) * 100, 2)) + PvalMiranda[species][domain]]
                               
        newfile.write('\t'.join(line) + '\n')        
newfile.close()
