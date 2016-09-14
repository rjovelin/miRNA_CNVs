# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 15:17:19 2016

@author: RJovelin
"""


# use this script to generate a table with target sites in CNV and non-CNv genes in each DGV release

# usage MakeTableTargetSitesDGCRelease.py [options]
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

CNVFiles = ['H_sapiens_GRCh37_2013-05_CNV_all_length_valid_chromos.txt',
            'H_sapiens_GRCh37_2013-07_CNV_all_length_valid_chromos.txt',
            'H_sapiens_GRCh37_2014_CNV_all_length_valid_chromos.txt',
            'H_sapiens_GRCh37_2015_CNV_all_length_valid_chromos.txt'] 

# make a dictionnary {release: {gene: [N_targets, Sequence_length, N_targets_normalized, CNV_status}}
# for each predictor, targetscan and miranda
targetscan, miranda = {}, {}
predictors = ['targetscan', 'miranda']

# get the seq input file
seq_input_file = 'H_sapiens_' + domain + '_' + chromos + '_targetscan.txt'
print(seq_input_file)

# loop over predictors
for i in range(len(predictors)):
    # get the predicted targets output file
    predicted_targets = 'H_sapiens_' + domain + '_' + chromos + '_predicted_sites_' + predictors[i] + '.txt'
    print(predicted_targets)
    # parse the predictor outputfile to get a dict {gene: [targets, seq_length, normalized_targets]}
    if i == 0:
        PredictedTargets = parse_targetscan_output(seq_input_file, predicted_targets, 'all')
    elif i == 1:
        PredictedTargets = parse_miranda_output(seq_input_file, predicted_targets, 'all')
    print('targets', len(PredictedTargets))    
    
    # loop over CNVFiles
    for filename in CNVFiles:
        print(filename)
        # get the release version
        release = filename[filename.index('GRCh37_')+len('GRCh37_'): filename.index('_CNV_')]
        if release == '2013_05':
            release = '2013a'
        elif release == '2013_07':
            release = '2013b'
        
        # get CNV gene status
        CNV_status = sort_genes_CNV_status(filename)
        print('CNV status', len(CNV_status))
        
        # copy dictionary with predicted targets
        targets = dict(PredictedTargets)
               
        # add CNV status
        for gene in targets:
            targets[gene].append(CNV_status[gene])
            
        # initialize inner dict
        if i == 0:
            targetscan[release] = {}
            for gene in targets:
                targetscan[release][gene] = list(targets[gene])
            print(release, len(targetscan[release][list(targetscan[release].keys())[0]]))
            print('targets', len(PredictedTargets[list(PredictedTargets.keys())[0]]))
        elif i == 1:
            miranda[release] = {}
            for gene in targets:
                miranda[release][gene] = list(targets[gene])
            print(release, len(miranda[release][list(miranda[release].keys())[0]]))
            print('targets', len(PredictedTargets[list(PredictedTargets.keys())[0]]))
        
# check that the same number of species are recorded for miranda and targetscan
assert len(targetscan) == len(miranda), 'different number of species depending on predictor'
# check that list values have the correct number of items
for release in targetscan:
    for gene in targetscan[release]:
        assert len(targetscan[release][gene]) == 4, 'gene in targetscan does not have all required values'
    for gene  in miranda[release]:
        assert len(miranda[release][gene]) == 4, 'gene in miranda does not have all required values'


# create a dicts with species as keys and list of targets for CNV and non-CNV genes,
# {release : [[CNV], [non-CNV]]}
DataTargetscan, DataMiranda  = {}, {}

# loop over release, get the number of sites for CNV and non-CNV genes
for release in targetscan:
    # initialise list value
    DataTargetscan[release] = [[], []]
    # populate inner lists with number of miRNA target sites
    for gene in targetscan[release]:
        if counts == 'raw' and targetscan[release][gene][-1] == 'CNV':
            DataTargetscan[release][0].append(targetscan[release][gene][0])
        elif counts == 'normalized' and targetscan[release][gene][-1] == 'CNV':
            DataTargetscan[release][0].append(targetscan[release][gene][2])
        elif counts == 'raw' and targetscan[release][gene][-1] == 'not_CNV':
            DataTargetscan[release][1].append(targetscan[release][gene][0])
        elif counts == 'normalized' and targetscan[release][gene][-1] == 'not_CNV':
            DataTargetscan[release][1].append(targetscan[release][gene][2])
for species in miranda:
    # initialize list values
    DataMiranda[release] = [[], []]
    # populate inner lists with number of mirna target sites
    for gene in miranda[release]:
        if counts == 'raw' and miranda[release][gene][-1] == 'CNV':
            DataMiranda[release][0].append(miranda[release][gene][0])
        elif counts == 'normalized' and miranda[release][gene][-1] == 'CNV':
            DataMiranda[release][0].append(miranda[release][gene][2])
        elif counts == 'raw' and miranda[release][gene][-1] == 'not_CNV':
            DataMiranda[release][1].append(miranda[release][gene][0])
        elif counts == 'normalized' and miranda[release][gene][-1] == 'not_CNV':
            DataMiranda[release][1].append(miranda[release][gene][2])
print('generated lists of target sites for CNV and non-CNV genes')

# generate dicts with gene counts {release: [N_CNV_genes, N_nonCNV_genes]}
GeneNumbers = {}
for release in targetscan:
    # initialise counters
    i, j, k, l = 0, 0, 0, 0
    # update counters
    for gene in targetscan[release]:
        if targetscan[release][gene][-1] == 'CNV':
            i += 1
        elif targetscan[release][gene][-1] == 'not_CNV':
            j += 1
    for gene in miranda[release]:
        if miranda[release][gene][-1] == 'CNV':
            k += 1
        elif miranda[release][gene][-1] == 'not_CNV':
            l += 1
    # check that numbers match
    assert i == k, 'CNV genes should match between miranda and targetscan'
    assert j == l, 'non-CNV genes should match between miranda and targetscan'
    # populate dict
    GeneNumbers[release] = [i, j]
print('counted CNV and non-CNV genes for each species')

# perform statistical tests between CNV and non-CNV genes
# create dicts to store results {release: [P-value targetscan, P-value mirnada]}
CompTargets = {}
for release in DataTargetscan:
    Ptargetscan = stats.ranksums(DataTargetscan[release][0], DataTargetscan[release][1])[1]
    Pmiranda = stats.ranksums(DataMiranda[release][0], DataMiranda[release][1])[1]
    CompTargets[release] = [Ptargetscan, Pmiranda]
print('compared CNV and non-CNV genes')
# get the significance level
Significance = {}
for release in CompTargets:
    Significance[release] = []
    for i in range(len(CompTargets[release])):
        if CompTargets[release][i] >= 0.05:
            Significance[release].append('')
        elif CompTargets[release][i] < 0.05 and CompTargets[release][i] >= 0.01:
            Significance[release].append('*')
        elif CompTargets[release][i] < 0.01 and CompTargets[release][i] >= 0.001:
            Significance[release].append('**')
        elif CompTargets[release][i] < 0.001:
            Significance[release].append('***')


# write results to file
# get outputfile name
if counts == 'raw':
    outputfile = 'TableDGVRawCounts_' + domain + '_' + chromos + '_' + cnv_length + '.txt'
elif counts == 'normalized':
    outputfile = 'TableDGVNormalizedCounts_' + domain + '_' + chromos + '_' + cnv_length + '.txt'
print(outputfile)

# open file for writing
newfile = open(outputfile, 'w')
# write table header
newfile.write('\t'.join(['', '', '', 'TargetScan', '', '', 'miRanda', '', '']) + '\n')
newfile.write('\t'.join(['', 'CNV', 'Non-CNV', 'CNV', 'Non-CNV', '', 'CNV', 'Non-CNV', '']) + '\n')
newfile.write('\t'.join(['Sp', 'Na', 'Na', 'Mean', 'Mean', 'D (%)b', 'Mean', 'Mean', 'D (%)b']) + '\n')

# make a list of DGV releases loop over
Versions = ['2013a', '2013b', '2014', '2015']

for release in Versions:
    # create the line to write 
    line = [release, str(GeneNumbers[release][0]), str(GeneNumbers[release][1]),
            str(round(np.mean(DataTargetscan[release][0]), 4)), str(round(np.mean(DataTargetscan[release][1]), 4)),
            str(round((1 - np.mean(DataTargetscan[release][1])/np.mean(DataTargetscan[release][0])) * 100, 2)) + Significance[release][0],
            str(round(np.mean(DataMiranda[release][0]), 4)), str(round(np.mean(DataMiranda[release][1]), 4)),
            str(round((1 - np.mean(DataMiranda[release][1])/np.mean(DataMiranda[release][0])) * 100, 2)) + Significance[release][1]] 
    
    newfile.write('\t'.join(line) + '\n')

newfile.close()
