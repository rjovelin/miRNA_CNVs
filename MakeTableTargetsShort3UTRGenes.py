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

# loop over predictors
for i in range(len(predictors)):
    # loop over species names
    for species in species_codes:
        # loop over domain
        for domain in regions:
            print(predictors[i], species, domain)      
        
        
        
        
        
        
        
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
            targetscan[species] = {}
            for gene in targets:
                # record gene with short 3'UTR
                if gene in UTR_length and UTR_length[gene] == 'short':
                    targetscan[species][gene] = list(targets[gene])
        elif i == 1:
            miranda[species] = {}
            for gene in targets:
                # record gene with short 3'UTR
                if gene in UTR_length and UTR_length[gene] == 'short':
                    miranda[species][gene] = list(targets[gene])



        
# check that the same number of species are recorded for miranda and targetscan
assert len(targetscan) == len(miranda), 'different number of species depending on predictor'
# check that list values have the correct number of items
for species in targetscan:
    for gene in targetscan[species]:
        assert len(targetscan[species][gene]) == 4, 'gene in targetscan does not have all required values'
    for gene  in miranda[species]:
        assert len(miranda[species][gene]) == 4, 'gene in miranda does not have all required values'


# create a dict of {species name : [[CNV], [non-CNV]]}
SpeciesDataTargetscan, SpeciesDataMiranda = {}, {}

# loop over species names, get the number of normalized sites for CNV and non-CNV genes
for species in targetscan:
    # initialise list value
    SpeciesDataTargetscan[species] = [[], []]
    # populate inner lists with number of miRNA target sites per nucleotide
    for gene in targetscan[species]:
        if targetscan[species][gene][-1] == 'CNV':
            SpeciesDataTargetscan[species][0].append(targetscan[species][gene][2])
        elif targetscan[species][gene][-1] == 'not_CNV':
            SpeciesDataTargetscan[species][1].append(targetscan[species][gene][2])
for species in miranda:
    # initialize list values
    SpeciesDataMiranda[species] = [[], []]
    # populate inner lists with number of mirna target sites per nucleotide
    for gene in miranda[species]:
        if miranda[species][gene][-1] == 'CNV':
            SpeciesDataMiranda[species][0].append(miranda[species][gene][2])
        elif miranda[species][gene][-1] == 'not_CNV':
            SpeciesDataMiranda[species][1].append(miranda[species][gene][2])
print('generated lists of target sites for CNV and non-CNV genes')


# perform stattistical tests between CNV and non-CNV genes
# create a dict to store results {species: P-value}
CompTargetscan, CompMiranda = {}, {}
for species in SpeciesDataTargetscan:
    P = stats.ranksums(SpeciesDataTargetscan[species][0], SpeciesDataTargetscan[species][1])[1]
    CompTargetscan[species] = P
for species in SpeciesDataMiranda:
    P = stats.ranksums(SpeciesDataMiranda[species][0], SpeciesDataMiranda[species][1])[1]    
    CompMiranda[species] = P
print('compared CNV and non-CNV genes')

# print P-values
for species in CompTargetscan:
    print('targetscan', species, CompTargetscan[species])
for species in CompMiranda:
    print('miranda', species, CompMiranda[species])


# make a list of data for each predictor
AllDataTargetscan, AllDataMiranda = [], []

# make a list of species names to loop from
species_names = ['H_sapiens', 'P_troglodytes', 'M_mulatta', 'M_musculus', 'B_taurus', 'G_gallus']
# loop over species in species names list and populate data lists, keeping the same order for targetscan and miranda
for species in species_names:
    # append list of target sites for CNV genes
    AllDataTargetscan.append(SpeciesDataTargetscan[species][0])
    AllDataMiranda.append(SpeciesDataMiranda[species][0])
    # append list of target sites for non-CNV genes
    AllDataTargetscan.append(SpeciesDataTargetscan[species][1])
    AllDataMiranda.append(SpeciesDataMiranda[species][1])
print('data consolidated in array')


# create figure
fig = plt.figure(1, figsize = (8, 3))

# create list of labels and tick positions for the X axis
#xtickpos = [0.35, 1.25, 2.15, 3.05, 3.95, 4.85]
xtickpos = [0.2, 1.1, 2, 2.9, 3.8, 4.7]
Names = [species_codes[i] for i in species_names]
print(Names)


# annotate Graph with significance level
PvalTargetScan, PvalMiranda = [], []
for species in species_names:
    if CompTargetscan[species] >= 0.05:
        PvalTargetScan.append('')
    elif CompTargetscan[species] < 0.05 and CompTargetscan[species] >= 0.01:
        PvalTargetScan.append('*')
    elif CompTargetscan[species] < 0.01 and CompTargetscan[species] >= 0.001:
        PvalTargetScan.append('**')
    elif CompTargetscan[species] < 0.001:
        PvalTargetScan.append('***')
for species in species_names:
    if CompMiranda[species] >= 0.05:
        PvalMiranda.append('')
    elif CompMiranda[species] < 0.05 and CompMiranda[species] >= 0.01:
        PvalMiranda.append('*')
    elif CompMiranda[species] < 0.01 and CompMiranda[species] >= 0.001:
        PvalMiranda.append('**')
    elif CompMiranda[species] < 0.001:
        PvalMiranda.append('***')



########################



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
            str(round((1 - np.mean(SpeciesDataMiranda[species][1])/np.mean(SpeciesDataMiranda[species][0])) * 100, 2)) + Significance[species][1]] 
    
    newfile.write('\t'.join(line) + '\n')

newfile.close()
