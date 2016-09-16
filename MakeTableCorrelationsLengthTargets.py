# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 15:48:57 2016

@author: RJovelin
"""


# use this script to generate a table with correlation coefficient between
# sequence length and number of targets for each species

# usage MakeTableCorrelationsLengthTargets.py

# import modules
import numpy as np
from scipy import stats
import math
import os
import sys
# import custom modules
from CNV_miRNAs import *


# keep genes on assembled nuclear chromosomes
chromos = 'valid_chromos'
keep_valid_chromos = True
# consider all CNVs
cnv_length = 'CNV_all_length'
CNV_size = 'all'
# alternatives
# chromos = 'all_chromos'
# cnv_length = 'CNV_greater_1Kb'
# CNV_size = 'long'

# make a dictionary of species names : species code
species_codes = {'H_sapiens': 'Hsa', 'P_troglodytes': 'Ptr', 'M_mulatta': 'Mml',
                 'M_musculus': 'Mmu', 'B_taurus': 'Bta', 'G_gallus':'Gga'}

regions = ['3UTR', '5UTR', 'CDS']
predictors = ['targetscan', 'miranda']
# make a dict {species: {domain: correlation}}
CorrelTargetscan, CorrelMiranda = {}, {}

# loop over species
for species in species_codes:
    # loop over domain
    for domain in regions:
        # loop over predictors
        for i in range(len(predictors)):
            print(species, domain, predictors[i])
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
            # make parallel lists with domain length and number of targets
            Length, Targets = [], []
            for gene in targets:
                Length.append(targets[gene][1])
                Targets.append(targets[gene][0])
            # take the spearman's rho correlation
            rho = stats.spearmanr(Length, Targets)[0]
            # check if predictor is targetscan or miranda
            if i == 0:
                # check if species in dict
                if species not in CorrelTargetscan:
                    CorrelTargetscan[species] = {}
                else:
                    CorrelTargetscan[species][domain] = rho
            elif i == 1:
                # check if species in dict
                if species not in CorrelMiranda:
                    CorrelMiranda[species] = {}
                else:
                    CorrelMiranda[species][domain] = rho
            
 
# write results to file
outputfile = 'TableCorrelationLengthTargets_' + domain + '_' + chromos + '_' + cnv_length + '.txt'
print(outputfile)

# make a list of species names to loop over
SpeciesNames = ['H_sapiens', 'P_troglodytes', 'M_mulatta', 'M_musculus', 'B_taurus', 'G_gallus']

# open file for writing
newfile = open(outputfile, 'w')
# write table header
newfile.write('\t'.join(['', 'TargetScan', '', '', 'miRanda', '', '']) + '\n')
newfile.write('\t'.join(['Species', '3\'UTR', '5\'UTR', 'CDS', '3\'UTR', '5\'UTR', 'CDS']) + '\n')

for species in SpeciesNames:
    # create the line to write to file
    line = [species_code[species]]
    for domain in regions:
        line.append(str(round(CorrelTargetscan[species][domain], 4)))
    for domain in regions:
        line.append(str(round(CorrelMiranda[species][domain], 4)))
    newfile.write('\t'.join(line) + '\n')         
newfile.close()
