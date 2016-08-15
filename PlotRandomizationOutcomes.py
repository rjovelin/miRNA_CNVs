# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 11:42:01 2016

@author: RJovelin
"""



# use this script to determine the proportion of random datasets with CNV genes having more, less or similar
# miRNA target sites as non-CNV genes

# usage PlotBootstrapOutcomes.py [options]
# [3UTR/5UTR/CDS]: gene domain to consider


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
import random
# import custom modules
from CNV_miRNAs import *


# get the region to consider to predict target sites [3UTR or 5UTr or CDS]
domain = sys.argv[1]
print(domain)

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
        
        # initialize inner dict
        if i == 0:
            targetscan[species] = {}
            for gene in targets:
                targetscan[species][gene] = list(targets[gene])
        elif i == 1:
            miranda[species] = {}
            for gene in targets:
                miranda[species][gene] = list(targets[gene])

print('obtained target sites')
# check that the same number of species are recorded for miranda and targetscan
assert len(targetscan) == len(miranda), 'different number of species depending on predictor'

# write functions used to randonly bootstrap genes

# use this function to match numbers and genes to randomly sample genes
def AddNumberToGenes(targets):
    '''
    (dict) -> dict
    Take the dict with target sites per gene and species and return a dict with
    number : gene pairs for each species    
    '''
    # create a dictionary {species: {num: gene}}} 
    tosamplefrom = {}
    # loop over species
    for species in targets:
        # initialize dict
        tosamplefrom[species] = {}
        # initialize counter
        i = 0
        # loop over genes for given species
        for gene in targets[species]:
            # populate dict with num: gene pair and update counter
            tosamplefrom[species][i] = gene
            i += 1
    return tosamplefrom


# create a function to randomzie gene CNV status
def RandomizeCNVStatus(tosamplefrom, targets):
    '''
    (dict, dict) -> dict
    Take a dictionary with number : gene pairs and a dictionary with mirna targets
    and return a dictionary with species as key and a list with number of random
    datasets for which CNV genes have more targets, less targets or similar
    number of targets as non-CNV genes
    Precondition: use normalized number of targets
    '''
    
    # tosamplefrom  is the form {species: {num: gene}}} 
    # targets is the form {species: {gene: [targets, seq_length, normalized_targets]}}    
    
    # create a dict for each species with a list with numbers of each different outcomes
    # when comparing targets in CNV and non-CNV genes
    # {species: [# replicates CNV > non-CNV, # replicates CNV < non-CNV, # replicates no differences]}
    RandomizedDatasets = {}
    # initialize list values
    for species in tosamplefrom:
        RandomizedDatasets[species] = [0, 0, 0]
    # create a list with CNV status to randomly assign status to genes
    CNVStatus = ['CNV', 'not_CNV']
    # loop over studies in dict to sample from
    for species in tosamplefrom:
        print('randomizing', species)
        # set number of replicates
        replicates = 10000
        # create lists with number of cnv and non-cnv genes for each random dataset
        a, b = [], []
        while replicates != 0:
            # make list of targets for CNV and non-CNV genes
            repCNVtargets, repNonCNVtargets = [], []
            # draw 1000 genes at random
            for i in range(1000):
                # draw a random CNV gene
                j = random.randint(0, len(tosamplefrom[species]) - 1)
                # get the corresponding genes
                randomgene = tosamplefrom[species][j]
                # assign CNV status at random
                k = random.randint(0, 1)
                RandomCNVStatus = CNVStatus[k]
                # check if gene is CNV or not CNV
                if RandomCNVStatus == 'CNV':
                    repCNVtargets.append(targets[species][randomgene][2])
                elif RandomCNVStatus == 'not_CNV':
                    repNonCNVtargets.append(targets[species][randomgene][2])
            # make sure that the correct numbers of genes is drawn
            assert len(repCNVtargets) + len(repNonCNVtargets) == 1000, '1000 non-CNV genes should be drawn'
            # record the number of cnv and non-cnv genes            
            a.append(len(repCNVtargets))
            b.append(len(repNonCNVtargets))
            # compare CNV and non-CNV genes
            Pval = stats.ranksums(repCNVtargets, repNonCNVtargets)[1]
            # check significance
            if Pval >= 0.05:
                # difference is not significant
                RandomizedDatasets[species][2] += 1
            elif Pval < 0.05:
                # difference is significance, check if CNV genes have a greater number of targets
                if np.mean(repCNVtargets) > np.mean(repNonCNVtargets):
                    RandomizedDatasets[species][0] += 1
                elif np.mean(repCNVtargets) < np.mean(repNonCNVtargets):
                    RandomizedDatasets[species][1] += 1
            # update replicate number
            replicates -= 1
        # print mean number of cnv and non-cnv genes
        print(species, 'mean CNV genes', np.mean(a), 'mean non-cnv genes', np.mean(b))

    return RandomizedDatasets


# assign numbers to gene to sample from {species: {num: gene}}} 
ToSampleFromTargetScan = AddNumberToGenes(targetscan)
ToSampleFromMiranda = AddNumberToGenes(miranda)
print('matched genes with numbers for sampling')        

RandomizedTargetscan = RandomizeCNVStatus(ToSampleFromTargetScan, targetscan)
print('done randomizing for targetscan')        
RandomizedMiranda = RandomizeCNVStatus(ToSampleFromMiranda, miranda)
print('done randomizing for miranda')        

# print outcomes
for species in SpeciesNames:
    print(RandomizedTargetscan[species])
for species in SpeciesNames:
    print(RandomizedMiranda[species])

# create parallel list of proportions 
GreaterTargetscan, LowerTargetscan, NodiffTargetscan, GreaterMiranda, LowerMiranda, NodiffMiranda = [], [] ,[], [], [], []
SpeciesNames = ['H_sapiens', 'P_troglodytes', 'M_mulatta', 'M_musculus', 'B_taurus', 'G_gallus']
# loop over species, get the proportions of replicates for each outcome
for species in SpeciesNames:
    GreaterTargetscan.append(RandomizedTargetscan[species][0] / sum(RandomizedTargetscan[species]))
    LowerTargetscan.append(RandomizedTargetscan[species][1] / sum(RandomizedTargetscan[species]))
    NodiffTargetscan.append(RandomizedTargetscan[species][2] / sum(RandomizedTargetscan[species]))
    GreaterMiranda.append(RandomizedMiranda[species][0] / sum(RandomizedMiranda[species]))
    LowerMiranda.append(RandomizedMiranda[species][1] / sum(RandomizedMiranda[species]))
    NodiffMiranda.append(RandomizedMiranda[species][2] / sum(RandomizedMiranda[species]))

# create a list with all the proportion lists
ProportionsTargetscan = [GreaterTargetscan, LowerTargetscan, NodiffTargetscan]
ProportionsMiranda = [GreaterMiranda, LowerMiranda, NodiffMiranda]

# create figure
fig = plt.figure(1, figsize = (4, 6))

# create a function to format the subplots
def CreateAx(Columns, Rows, Position, Data, figure, Title, LabelNames, XScale):
    '''
    (int, int, int, list, figure_object, str, int, list, list)
    Take the number of a column, and rows in the figure object and the position of
    the ax in figure, a list of data, a title, a maximum value for the Y axis,
    a list with species names and list of X axis tick positions and return an
    ax instance in the figure
    '''    
    
    # create subplot in figure
    # add a plot to figure (N row, N column, plot N)
    ax = figure.add_subplot(Rows, Columns, Position)
    
    # get the list of proportions 
    greater, lower, nodiff = Data[0], Data[1], Data[2]
    # make a list with added values for nodiff and greater
    added = []
    for i in range(len(greater)):
        added.append(nodiff[i] + greater[i])
    # Create a bar plot for proportions of replicates with CNV no diff on top of CNV lower
    ax.bar([0, 0.4, 0.8, 1.2, 1.6, 2], nodiff, width = 0.3, label = 'No difference', color= '#f7f7f7')
    # Create a bar plot for proportions of replicates with CNV greater on top of no diff
    ax.bar([0, 0.4, 0.8, 1.2, 1.6, 2], greater, width = 0.3, bottom = nodiff, label = 'CNV > non-CNV', color= '#ef8a62')
    # Create a bar plot for proportions of replicates with CNV lower on top of CNV greater
    ax.bar([0, 0.4, 0.8, 1.2, 1.6, 2], lower, width = 0.3, bottom= added, label = 'CNV < non-CNV', color = '#67a9cf')
 
    # write title
    if Title == 'TargetScan':
        ax.set_title(Title + '\n\n', size = 10)
    elif Title == 'miRanda':
        ax.set_title(Title, size = 10)
    # set font for all text in figure
    FigFont = {'fontname':'Arial'}   
    # write label for y axis
    ax.set_ylabel('Proportion of replicates', color = 'black', size = 10, ha = 'center', **FigFont)
    # write label for x axis
    plt.xticks(XScale, LabelNames, ha = 'center', fontsize = 10, **FigFont)
        
    # do not show lines around figure  
    ax.spines["top"].set_visible(False)    
    ax.spines["bottom"].set_visible(True)    
    ax.spines["right"].set_visible(False)    
    ax.spines["left"].set_visible(False)  
    # offset the spines
    for spine in ax.spines.values():
        spine.set_position(('outward', 5))
    
    # add a light grey horizontal grid to the plot, semi-transparent, 
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5, linewidth = 0.5)  
    # hide these grids behind plot objects
    ax.set_axisbelow(True)

    # do not show ticks
    plt.tick_params(
        axis='both',       # changes apply to the x-axis and y-axis (other option : x, y)
        which='both',      # both major and minor ticks are affected
        bottom='on',      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        right = 'off',
        left = 'off',          
        labelbottom='on', # labels along the bottom edge are on
        colors = 'black',
        labelsize = 10,
        direction = 'out') # ticks are outside the frame when bottom = 'on'  
      
    # Set the tick labels font name
    for label in ax.get_yticklabels():
        label.set_fontname('Arial')
    
    # create a margin around the x axis
    plt.margins(0.05)
    
    return ax      

# plot porportions for targetscan and miranda
Names = [species_codes[i] for i in SpeciesNames]
ax1 = CreateAx(1, 2, 1, ProportionsTargetscan, fig, 'TargetScan', Names, [0.15, 0.55, 0.95, 1.35, 1.75, 2.15])
ax2 = CreateAx(1, 2, 2, ProportionsMiranda, fig, 'miRanda', Names, [0.15, 0.55, 0.95, 1.35, 1.75, 2.15])

# add subplot label
ax1.text(-0.5, 1.1, 'A', horizontalalignment = 'center',
         verticalalignment = 'center', color = 'black', size = 10)
ax2.text(-0.5, 1.1, 'B', horizontalalignment = 'center',
         verticalalignment = 'center', color = 'black', size = 10)

# add legend
N = mpatches.Patch(facecolor = '#f7f7f7' , edgecolor = 'black', linewidth = 1, label= 'No diff.')
G = mpatches.Patch(facecolor = '#ef8a62' , edgecolor = 'black', linewidth = 1, label= 'CNV greater')
L = mpatches.Patch(facecolor = '#67a9cf' , edgecolor = 'black', linewidth = 1, label= 'CNV lower')
ax1.legend(handles = [N, G, L], loc = (0, 1), fontsize = 8, frameon = False, ncol = 3)

# make sure subplots do not overlap
plt.tight_layout()

# get outputfile
outputfile = 'PlotRandomization' + '_' + domain + '_' + chromos + '_' + cnv_length 
print(outputfile)

# save figure
fig.savefig(outputfile + '.eps', bbox_inches = 'tight')
       
