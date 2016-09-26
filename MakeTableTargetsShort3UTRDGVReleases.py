# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 11:17:34 2016

@author: RJovelin
"""



# use this script to generate a table with target sites in CNV and non-CNv genes
# with short 3'UTRs in each DGV release

# usage MakeTableTargetShort3UTRDGVReleases.py [options]
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

# make a dictionnary {release: {domain: {gene: [N_targets, Sequence_length, N_targets_normalized, CNV_status}}}
# for each predictor, targetscan and miranda
targetscan, miranda = {}, {}
predictors = ['targetscan', 'miranda']
regions = ['5UTR', 'CDS']


# loop over CNV files
for filename in CNVFiles:
    # loop over predictors
    for i in range(len(predictors)):
        # loop over domain
        for domain in regions:
            print(i, predictors[i], filename, domain)
            # get the seq input file
            seq_input_file = 'H_sapiens_' + domain + '_' + chromos + '_targetscan.txt'
            print(seq_input_file)
            # get the predicted targets output file
            predicted_targets = 'H_sapiens_' + domain + '_' + chromos + '_predicted_sites_' + predictors[i] + '.txt'
            print(predicted_targets)
            # parse the predictor outputfile to get a dict {gene: [targets, seq_length, normalized_targets]}
            if i == 0:
                targets = parse_targetscan_output(seq_input_file, predicted_targets, 'all')
            elif i == 1:
                targets = parse_miranda_output(seq_input_file, predicted_targets, 'all')
            print('targets', len(targets))    

            # get the release version
            release = filename[filename.index('GRCh37_')+len('GRCh37_'): filename.index('_CNV_')]
            if release == '2013-05':
                release = '2013a'
            elif release == '2013-07':
                release = '2013b'
            # get CNV gene status
            CNV_status = sort_genes_CNV_status(filename)
            print('CNV status', len(CNV_status))
            # add CNV status
            for gene in targets:
                targets[gene].append(CNV_status[gene])
            
            # get file with UTR status (short, long)
            UTR_file = 'H_sapiens_3UTR_length_' + chromos + '.txt'
            print(UTR_file)
            # create a dict {gene : "short" (or "long")}
            UTR_length = sort_genes_3UTR_length(UTR_file, L)
            print('UTR length', len(UTR_length))
            # initialize inner dict
            if i == 0:
                # check if species in dict
                if release not in targetscan:
                    targetscan[release] = {}
                #initialize inner dict with domain as key
                targetscan[release][domain] = {}
                for gene in targets:
                    # record gene with short 3'UTR
                    if gene in UTR_length and UTR_length[gene] == 'short':
                         targetscan[release][domain][gene] = copy.deepcopy(targets[gene])
            elif i == 1:
                # check if species in dict
                if release not in miranda:
                    miranda[release] = {}
                # initialize inner dict with domain as key
                miranda[release][domain] = {}
                for gene in targets:
                    # record gene with short 3'UTR
                    if gene in UTR_length and UTR_length[gene] == 'short':
                        miranda[release][domain][gene] = copy.deepcopy(targets[gene])
            
# check that the same number of species are recorded for miranda and targetscan
assert len(targetscan) == len(miranda), 'different number of species depending on predictor'
# check that list values have the correct number of items
for release in targetscan:
    for domain in targetscan[release]:
        for gene in targetscan[release][domain]:
            assert len(targetscan[release][domain][gene]) == 4, 'gene in targetscan does not have all required values'
        for gene  in miranda[release][domain]:
            assert len(miranda[release][domain][gene]) == 4, 'gene in miranda does not have all required values'


# create a dict of {species name : {domain: [CNV], [non-CNV]]}}
DataTargetscan, DataMiranda = {}, {}

# loop over DGV releases, get the number of normalized sites for CNV and non-CNV genes
for release in targetscan:
    # initialise inner dict
    DataTargetscan[release] = {}
    # loop over domain
    for domain in targetscan[release]:
        # initialize list value
        DataTargetscan[release][domain] = [[], []]    
        # populate inner lists with number of miRNA target sites per nucleotide
        for gene in targetscan[release][domain]:
            if targetscan[release][domain][gene][-1] == 'CNV':
                DataTargetscan[release][domain][0].append(targetscan[release][domain][gene][2])
            elif targetscan[release][domain][gene][-1] == 'not_CNV':
                DataTargetscan[release][domain][1].append(targetscan[release][domain][gene][2])
for release in miranda:
    # initialize inner dict
    DataMiranda[release] = {}
    # loop over domain
    for domain in miranda[release]:
        # initialize list value
        DataMiranda[release][domain] = [[], []]
        # populate inner lists with number of mirna targets per nucleotide
        for gene in miranda[release][domain]:
            if miranda[release][domain][gene][-1] == 'CNV':
                DataMiranda[release][domain][0].append(miranda[release][domain][gene][2])
            elif miranda[release][domain][gene][-1] == 'not_CNV':
                DataMiranda[release][domain][1].append(miranda[release][domain][gene][2])
print('generated lists of target sites for CNV and non-CNV genes')

# create a dict with gene numbers {release: {domain: [N CNV genes targetscan, N CNV genes miranda]}}
GeneNumbersTargetscan, GeneNumbersMiranda = {}, {}
for release in targetscan:
    GeneNumbersTargetscan[release] = {}
    for domain in targetscan[release]:
        GeneNumbersTargetscan[release][domain] = [0, 0]
        for gene in targetscan[release][domain]:
            if targetscan[release][domain][gene][-1] == 'CNV':
                GeneNumbersTargetscan[release][domain][0] += 1
            elif targetscan[release][domain][gene][-1] == 'not_CNV':
                GeneNumbersTargetscan[release][domain][1] += 1
for release in miranda:
    GeneNumbersMiranda[release] = {}
    for domain in miranda[release]:
        GeneNumbersMiranda[release][domain] = [0, 0]
        for gene in miranda[release][domain]:
            if miranda[release][domain][gene][-1] == 'CNV':
                GeneNumbersMiranda[release][domain][0] += 1
            elif miranda[release][domain][gene][-1] == 'not_CNV':
                GeneNumbersMiranda[release][domain][1] += 1
print('counted CNV and non-CNV genes')
            
# perform stattistical tests between CNV and non-CNV genes
# create a dict to store results {release: {domain: significance level}}
PvalTargetscan, PvalMiranda = {}, {}
for release in DataTargetscan:
    PvalTargetscan[release], PvalMiranda[release] = {}, {}
    for domain in DataTargetscan[release]:
        P = stats.ranksums(DataTargetscan[release][domain][0], DataTargetscan[release][domain][1])[1]
        if P >= 0.05:
            PvalTargetscan[release][domain] = ''
        elif P < 0.05 and P >= 0.01:
            PvalTargetscan[release][domain] = '*'
        elif P < 0.01 and P >= 0.001:
            PvalTargetscan[release][domain] = '**'
        elif P < 0.001:
            PvalTargetscan[release][domain] = '***'
    for domain in DataMiranda[release]:
        P = stats.ranksums(DataMiranda[release][domain][0], DataMiranda[release][domain][1])[1]
        if P >= 0.05:
            PvalMiranda[release][domain] = ''
        elif P < 0.05 and P >= 0.01:
            PvalMiranda[release][domain] = '*'
        elif P < 0.01 and P >= 0.001:
            PvalMiranda[release][domain] = '**'
        elif P < 0.001:
            PvalMiranda[release][domain] = '***'
print('compared CNV and non-CNV genes')


# write output to file
# get outputfile name
outputfile = 'TableTargetsShort3UTRDGV_' + str(L) + 'bp_' + chromos + '_' + cnv_length + '.txt'
# open file for writing
newfile = open(outputfile, 'w')
# write table header
newfile.write('\t'.join(['', '', '', '', 'TargetScan', '', '', 'miRanda', '', '']) + '\n')
newfile.write('\t'.join(['', '', 'CNV', 'Non-CNV', 'CNV', 'Non-CNV', '', 'CNV', 'Non-CNV', '']) + '\n')
newfile.write('\t'.join(['Release', 'Domain', 'Na', 'Na', 'Mean', 'Mean', 'D (%)b', 'Mean', 'Mean', 'D (%)b']) + '\n')


# make a list of DGV releases loop over
Versions = ['2013a', '2013b', '2014', '2015']

# loop over domains
for domain in regions:
    # loop over releases
    for release in Versions:
        # create the line to write
        line = [release, domain, str(GeneNumbersTargetscan[release][domain][0]), str(GeneNumbersTargetscan[release][domain][1]),
                str(round(np.mean(DataTargetscan[release][domain][0]), 4)), str(round(np.mean(DataTargetscan[release][domain][1]), 4)),
                str(round((1 - np.mean(DataTargetscan[release][domain][1])/np.mean(DataTargetscan[release][domain][0])) * 100, 2)) + PvalTargetscan[release][domain],
                str(round(np.mean(DataMiranda[release][domain][0]), 4)), str(round(np.mean(DataMiranda[release][domain][1]), 4)),
                str(round((1 - np.mean(DataMiranda[release][domain][1])/np.mean(DataMiranda[release][domain][0])) * 100, 2)) + PvalMiranda[release][domain]]
        newfile.write('\t'.join(line) + '\n')        
newfile.close()
