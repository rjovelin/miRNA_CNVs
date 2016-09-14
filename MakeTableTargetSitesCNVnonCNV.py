# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 13:46:09 2016

@author: RJovelin
"""


# use this script to generate a table with target sites in CNV and non-CNv genes in each species

# usage MakeTableTargetSitesCNVnonCNV.py [options]
# [3UTR/5UTR/CDS]: gene domain to consider
# [raw/normalized]: use raw or normalized target counts

# import modules
import numpy as np
from scipy import stats
import math
import os
import sys
# import custom modules
from CNV_miRNAs import *


# get the region to consider to predict target sites [3UTR or 5UTr or CDS]
domain = sys.argv[1]
# use raw target counts or normalized counts
counts = sys.argv[2]
print(domain, counts)

assert domain in ['3UTR', '5UTR', 'CDS'], 'domain is not recognized'
assert counts in ['raw', 'normalized'], 'use raw or normalized counts'


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

# make a dictionnary {species: {gene: [N_targets, Sequence_length, N_targets_normalized, CNV_status}}
# for each predictor, targetscan and miranda
targetscan, miranda = {}, {}

predictors = ['targetscan', 'miranda']
# loop over predictors
for i in range(len(predictors)):
    # loop over species names
    for species in species_codes:
        print(species)
        
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
            
        # initialize inner dict
        if i == 0:
            targetscan[species] = {}
            for gene in targets:
                targetscan[species][gene] = list(targets[gene])
        elif i == 1:
            miranda[species] = {}
            for gene in targets:
                miranda[species][gene] = list(targets[gene])
        
# check that the same number of species are recorded for miranda and targetscan
assert len(targetscan) == len(miranda), 'different number of species depending on predictor'
# check that list values have the correct number of items
for species in targetscan:
    for gene in targetscan[species]:
        assert len(targetscan[species][gene]) == 4, 'gene in targetscan does not have all required values'
    for gene  in miranda[species]:
        assert len(miranda[species][gene]) == 4, 'gene in miranda does not have all required values'


# create a dicts with species as keys and list of targets for CNV and non-CNV genes,
# {species name : [[CNV], [non-CNV]]}
SpeciesDataTargetscan, SpeciesDataMiranda  = {}, {}

# loop over species names, get the number of sites for CNV and non-CNV genes
for species in targetscan:
    # initialise list value
    SpeciesDataTargetscan[species] = [[], []]
    # populate inner lists with number of miRNA target sites
    for gene in targetscan[species]:
        if counts == 'raw' and targetscan[species][gene][-1] == 'CNV':
            SpeciesDataTargetscan[species][0].append(targetscan[species][gene][0])
        elif counts == 'normalized' and targetscan[species][gene][-1] == 'CNV':
            SpeciesDataTargetscan[species][0].append(targetscan[species][gene][2])
        elif counts == 'raw' and targetscan[species][gene][-1] == 'not_CNV':
            SpeciesDataTargetscan[species][1].append(targetscan[species][gene][0])
        elif counts == 'normalized' and targetscan[species][gene][-1] == 'not_CNV':
            SpeciesDataTargetscan[species][1].append(targetscan[species][gene][2])
for species in miranda:
    # initialize list values
    SpeciesDataMiranda[species] = [[], []]
    # populate inner lists with number of mirna target sites
    for gene in miranda[species]:
        if counts == 'raw' and miranda[species][gene][-1] == 'CNV':
            SpeciesDataMiranda[species][0].append(miranda[species][gene][0])
        elif counts == 'normalized' and miranda[species][gene][-1] == 'CNV':
            SpeciesDataMiranda[species][0].append(miranda[species][gene][2])
        elif counts == 'raw' and miranda[species][gene][-1] == 'not_CNV':
            SpeciesDataMiranda[species][1].append(miranda[species][gene][0])
        elif counts == 'normalized' and miranda[species][gene][-1] == 'not_CNV':
            SpeciesDataMiranda[species][1].append(miranda[species][gene][2])
print('generated lists of target sites for CNV and non-CNV genes')

# generate dicts with gene counts {species name: [N_CNV_genes, N_nonCNV_genes]}
GeneNumbers = {}
for species in targetscan:
    # initialise counters
    i, j, k, l = 0, 0, 0, 0
    # update counters
    for gene in targetscan[species]:
        if targetscan[species][gene][-1] == 'CNV':
            i += 1
        elif targetscan[species][gene][-1] == 'not_CNV':
            j += 1
    for gene in miranda[species]:
        if miranda[species][gene][-1] == 'CNV':
            k += 1
        elif miranda[species][gene][-1] == 'not_CNV':
            l += 1
    # check that numbers match
    assert i == k, 'CNV genes should match between miranda and targetscan'
    assert j == l, 'non-CNV genes should match between miranda and targetscan'
    # populate dict
    GeneNumbers[species] = [i, j]
print('counted CNV and non-CNV genes for each species')

# perform stattistical tests between CNV and non-CNV genes
# create dicts to store results {species: [P-value targetscan, P-value mirnada]}
CompTargets = {}
for species in SpeciesDataTargetscan:
    Ptargetscan = stats.ranksums(SpeciesDataTargetscan[species][0], SpeciesDataTargetscan[species][1])[1]
    Pmiranda = stats.ranksums(SpeciesDataMiranda[species][0], SpeciesDataMiranda[species][1])[1]
    CompTargets[species] = [Ptargetscan, Pmiranda]
print('compared CNV and non-CNV genes')
# get the significance level
Significance = {}
for species in CompTargets:
    Significance[species] = []
    for i in range(len(CompTargets[species])):
        if CompTargets[species][i] >= 0.05:
            Significance[species].append('')
        elif CompTargets[species][i] < 0.05 and CompTargets[species][i] >= 0.01:
            Significance[species].append('*')
        elif CompTargets[species][i] < 0.01 and CompTargets[species][i] >= 0.001:
            Significance[species].append('**')
        elif CompTargets[species][i] < 0.001:
            Significance[species].append('***')


# write results to file
# get outputfile name
if counts == 'raw':
    outputfile = 'TableRawCounts_' + domain + '_' + chromos + '_' + cnv_length + '.txt'
elif counts == 'normalized':
    outputfile = 'TableNormalizedCounts_' + domain + '_' + chromos + '_' + cnv_length + '.txt'
print(outputfile)

# open file for writing
newfile = open(outputfile, 'w')
# write table header
newfile.write('\t'.join(['', '', '', 'TargetScan', '', '', 'miRanda', '', '']) + '\n')
newfile.write('\t'.join(['', 'CNV', 'Non-CNV', 'CNV', 'Non-CNV', '', 'CNV', 'Non-CNV', '']) + '\n')
newfile.write('\t'.join(['Sp', 'Na', 'Na', 'Mean', 'Mean', 'D (%)b', 'Mean', 'Mean', 'D (%)b']) + '\n')

# make a list of species names to loop over
SpeciesNames = ['H_sapiens', 'P_troglodytes', 'M_mulatta', 'M_musculus', 'B_taurus', 'G_gallus']

for species in SpeciesNames:
    # create the line to write 
    line = [species_codes[species], str(GeneNumbers[species][0]), str(GeneNumbers[species][1]),
            str(round(np.mean(SpeciesDataTargetscan[species][0]), 4)), str(round(np.mean(SpeciesDataTargetscan[species][1]), 4)),
            str(round((1 - np.mean(SpeciesDataTargetscan[species][1])/np.mean(SpeciesDataTargetscan[species][0])) * 100, 2)) + Significance[species][0],
            str(round(np.mean(SpeciesDataMiranda[species][0]), 4)), str(round(np.mean(SpeciesDataMiranda[species][1]), 4)),
            str(round((1 - np.mean(SpeciesDataMiranda[species][1])/np.mean(SpeciesDataMiranda[species][0])) * 100), 2) + Significance[species][1]] 
    
    newfile.write('\t'.join(line) + '\n')

newfile.close()
