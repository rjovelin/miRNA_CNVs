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
    for domain in targetscan[species][domain]:
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
    for domain in miranda[species][domain]:
        DataMiranda[species][domain] = [[], []]
        # populate inner lists with number of mirna target sites per nucleotide
        for gene in miranda[species][domain]:
            if miranda[species][domain][gene][-1] == 'CNV':
                DataMiranda[species][domain][0].append(miranda[species][domain][gene][2])
            elif miranda[species][domain][gene][-1] == 'not_CNV':
                DataMiranda[species][domain][1].append(miranda[species][domain][gene][2])
print('generated lists of target sites for CNV and non-CNV genes')


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









Supplementary Table S8. Comparison of the mean number of miRNA binding sites predicted by TargetScan and miRanda normalized by the length of the 5’UTR sequences between CNV and non-CNV genes that have short 3’UTRs. aN: number of genes. bD: percent difference between CNV and non-CNV genes. A negative D indicates a greater mean number of target sites for non-CNV genes and a positive D indicates a greater mean number of target sites for CNV genes. * P < 0.05, ** P < 0.01, *** P < 0.001, Wilcoxon rank-sum tests. 
				TargetScan		miRanda
	CNV	Non-CNV		CNV	Non-CNV			CNV	Non-CNV	
Sp	Na	Na		Mean	Mean	D (%)b		Mean	Mean	D (%)b
Hsa	62	51		0.2567	0.2814	-8.78		0.1515	0.1733	-12.61
Ptr	17	733		0.0561	0.0576	-2.63		0.0337	0.0371	-9.28
Mml	264	1074		0.0863	0.0871	-0.90		0.0533	0.0540	-1.23
Mmu	51	285		0.1895	0.1939	-2.27		0.0968	0.1058	-8.52
Bta	68	577		0.0821	0.0819	0.22		0.0536	0.0467	14.78*
Gga	60	373		0.1113	0.1188	-6.27		0.0773	0.0751	2.88

Supplementary Table S9. Comparison of the mean number of miRNA binding sites predicted by TargetScan and miRanda normalized by the length of the CDS sequences between CNV and non-CNV genes that have short 3’UTRs. aN: number of genes. bD: percent difference between CNV and non-CNV genes. A negative D indicates a greater mean number of target sites for non-CNV genes and a positive D indicates a greater mean number of target sites for CNV genes. * P < 0.05, ** P < 0.01, *** P < 0.001, Wilcoxon rank-sum tests. 
				TargetScan		miRanda
	CNV	Non-CNV		CNV	Non-CNV			CNV	Non-CNV	
Sp	Na	Na		Mean	Mean	D (%)b		Mean	Mean	D (%)b
Hsa	473	140		0.2643	0.2958	-10.65***		0.1793	0.2041	-12.16***
Ptr	32	1813		0.0669	0.0669	-0.05		0.0462	0.0436	6.09
Mml	654	3760		0.1036	0.1010	2.54***		0.0633	0.0640	-1.05
Mmu	411	1443		0.2041	0.2054	-0.63		0.1289	0.1273	1.23
Bta	197	2151		0.0856	0.0897	-4.59***		0.0586	0.0587	-0.18
Gga	227	806		0.1317	0.1366	-3.59*		0.0837	0.0894	-6.37***


















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





########################


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
