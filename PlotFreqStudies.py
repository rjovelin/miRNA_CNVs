# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 16:53:07 2016

@author: RJovelin
"""


# use this script to plot the frequency of studies for which the number of
# miRNA sites for CNV genes is greater, lower or similar to non-CNV genes

# usage PlotFreqStudies.py [options]
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
# minimum_cnv =  minimum number of cnv genes in single study
MinimumCNVGenes = 500

Releases = ['GRCh37_2013-05', 'GRCh37_2013-07', 'GRCh37_2014', 'GRCh37_2015']
labelnames = ['2013a', '2013b', '2014', '2015']



# get the number of target sites for each study of each version of the DGV

# get UTR file
UTR_file = 'H_sapiens_3UTR_length_' + chromos + '.txt'
print(UTR_file)

# get synonym names for all genes {gene name : [list of synonyms]}
synonyms = get_synonyms('H_sapiens.gff3')
print('got synonymous names', len(synonyms))

# get the CDS sequences of the longest mRNAs for each gene {gene : sequence}
CDS_seq = extract_CDS_sequences('H_sapiens.gff3', 'H_sapiens_genome.txt', 'H_sapiens_valid_chromos.txt', keep_valid_chromos)
print('extracted CDS sequences', len(CDS_seq))  

# make a list of DGV files
DGVFiles = ['GRCh37_hg19_variants_2013-05-31.txt', 'GRCh37_hg19_variants_2013-07-23.txt',
            'GRCh37_hg19_variants_2014-10-16.txt', 'GRCh37_hg19_variants_2015-07-23.txt']

# make a dictionary to match each study to a pubmed ID in each release
# {release: {reference: pubmedid}}
References = {}
for filename in DGVFiles:
    # get release version
    if '2013' in filename:
        release_version = filename[:filename.index('_hg19')] + '_' + filename[filename.index('variants_') + len('variants_'): -7]
    else:
        release_version = filename[:filename.index('_hg19')] + '_' + filename[filename.index('variants_') + len('variants_'): -10]
    print(release_version)
    References[release_version] = {}
    ref = get_DGV_references(filename)
    References[release_version] = dict(ref)
print('got references')

# get the CNV genes for each study of each release of the DGV
# create a dict {release: {study: {set of cnv genes}}}
StudiesCNVGenes = {}
# loop over DGV files
for filename in DGVFiles:
    print(filename)
    # get release version
    if '2013' in filename:
        release_version = filename[:filename.index('_hg19')] + '_' + filename[filename.index('variants_') + len('variants_'): -7]
    else:
        release_version = filename[:filename.index('_hg19')] + '_' + filename[filename.index('variants_') + len('variants_'): -10]
    print(release_version)
    # initialize outer dict
    StudiesCNVGenes[release_version] = {}
    # get the set of CNV genes for each study of that release
    for study in References[release_version]:
        CNV_genes = get_human_CNV_genes_single_study(filename, study, 'all')
        StudiesCNVGenes[release_version][study] = set(CNV_genes)
print('got CNV genes for each study')        
        

# get the CNV status of all genes used to predict target sites for each study of each release
# create a dict {release: {study: {gene: CNV status}}}    
CNV_status = {}    
for release in StudiesCNVGenes:
    # initialize outer dict
    CNV_status[release] = {}
    for study in StudiesCNVGenes[release]:
        # initialize dict
        CNV_status[release][study] = {}
        # loop over genes in CDS seq
        for gene in CDS_seq:
            # set boolean
            is_cnv = False
            # ask if gene in CNV gene
            if gene in StudiesCNVGenes[release][study] or gene.upper() in StudiesCNVGenes[release][study]:
                CNV_status[release][study][gene] = 'CNV'
            else:
                # ask if any of the gene synonyms are in CNVs
                for name in synonyms[gene]:
                    # check if name in CNV
                    if name in StudiesCNVGenes[release][study] or name.upper() in StudiesCNVGenes[release][study]:
                        # update boolean
                        is_cnv = True
                # check if gene in CNV
                if is_cnv == True:
                    CNV_status[release][study][gene] = 'CNV'
                elif is_cnv == False:
                    CNV_status[release][study][gene] = 'not_CNV'
print('sorted genes according to CNV status')


# get the number of target sites for each gene
# get targetscan sequence input file
targetscan_seq_input_file = 'H_sapiens_' + domain + '_' + chromos + '_targetscan.txt'
print(targetscan_seq_input_file)
# get the outputfile with predicted target sites
TargetScanTargetsFile = 'H_sapiens_' + domain + '_' + chromos + '_predicted_sites_targetscan.txt'   
MirandaTargetsFile = 'H_sapiens_' + domain + '_' + chromos + '_predicted_sites_miranda.txt'   
# make a dictionary with {gene :[targets, seq_length, normalized_targets]}
TargetsTargetscan = parse_targetscan_output(targetscan_seq_input_file, TargetScanTargetsFile, 'all')
TargetsMiranda = parse_miranda_output(targetscan_seq_input_file, MirandaTargetsFile, 'all')

print('targetscan targets', len(TargetsTargetscan))
print('miranda targets', len(TargetsMiranda))

# create dicts for targetscan and miranda predictors with target sites and CNV status
# {release: {study: {gene: [targets, seq_length, normalized_targets, CNV_status]}}}
CNVTargetsTargetscan, CNVTargetsMiranda = {}, {}
# loop over release in CNV status
for release in CNV_status:
    # initialize inner dict
    CNVTargetsTargetscan[release], CNVTargetsMiranda[release] = {}, {}
    for study in CNV_status[release]:
        # initialize inner dict
        CNVTargetsTargetscan[release][study], CNVTargetsMiranda[release][study] = {}, {}
        # loop over genes, populate dict combining targets and CNv status 
        for gene in CNV_status[release][study]:
            # get CNV status and add it to the list of targets, without modifying the original list
            status = CNV_status[release][study][gene]
            if gene in TargetsTargetscan:
                CNVTargetsTargetscan[release][study][gene] = list(TargetsTargetscan[gene]) + [status]
            if gene in TargetsMiranda:
                CNVTargetsMiranda[release][study][gene] = list(TargetsMiranda[gene]) + [status]
print('combined CNV status and miRNA targets for each gene in each study')            


# compare the normalized number of target sites between CNV and non-CNV genes for each study and DGV version
# count the number of studies with CNV genes having greater, lowe or similar number targets as non-CNV genes
# {release: [CNV_greater, CNV_lower, NO_diff]}
CompTargetscan, CompMiranda = {}, {}
for release in CNVTargetsTargetscan:
    # initialize dict
    CompTargetscan[release] = [0, 0, 0]
    # loop over studies
    for study in CNVTargetsTargetscan[release]:
        # make lists of target sites for CNV amd non-CNV genes
        cnvtargets = [CNVTargetsTargetscan[release][study][gene][2] for gene in CNVTargetsTargetscan[release][study] if CNVTargetsTargetscan[release][study][gene][-1] == 'CNV']
        noncnvtargets = [CNVTargetsTargetscan[release][study][gene][2] for gene in CNVTargetsTargetscan[release][study] if CNVTargetsTargetscan[release][study][gene][-1] == 'not_CNV']
        # check that minimum number of cnv gene is met
        if len(cnvtargets) >= MinimumCNVGenes:
            # compare the number of target sites
            P = stats.ranksums(cnvtargets, noncnvtargets)[1]
            # update counters
            if P < 0.05:
                # compare means
                if np.mean(cnvtargets) > np.mean(noncnvtargets):
                    CompTargetscan[release][0] += 1
                elif np.mean(cnvtargets) < np.mean(noncnvtargets):
                    CompTargetscan[release][1] += 1
            elif P >= 0.05:
                CompTargetscan[release][2] += 1
print('compared mean targetscan sites between CNV and non-CNV genes')        
for release in CNVTargetsMiranda:
    # initialize dict
    CompMiranda[release] = [0, 0, 0]
    # loop over studies
    for study in CNVTargetsMiranda[release]:
        # make lists of target sites for CNV and non-CNV genes
        cnvtargets = [CNVTargetsMiranda[release][study][gene][2] for gene in CNVTargetsMiranda[release][study] if CNVTargetsMiranda[release][study][gene][-1] == 'CNV']
        noncnvtargets = [CNVTargetsMiranda[release][study][gene][2] for gene in CNVTargetsMiranda[release][study] if CNVTargetsMiranda[release][study][gene][-1] == 'not_CNV']
        # check that minimum number of cnv gene is met        
        if len(cnvtargets) >= MinimumCNVGenes:
            # compare the number of target sites
            P = stats.ranksums(cnvtargets, noncnvtargets)[1]
            # update counters
            if P < 0.05:
                # compare means
                if np.mean(cnvtargets) > np.mean(noncnvtargets):
                    CompMiranda[release][0] += 1
                elif np.mean(cnvtargets) < np.mean(noncnvtargets):
                    CompMiranda[release][1] += 1
            elif P >= 0.05:
                CompMiranda[release][2] += 1
print('compared mean miranda sites between CNV and non-CNV genes')        


# create parallel lists with CNV greater, lower and no diff for each release
CNVGreaterTargetscan, CNVLowerTargetscan, CNVNoDiffTargetscan = [], [], []
CNVGreaterMiranda, CNVLowerMiranda, CNVNoDiffMiranda = [], [], []
# populate lists with frequencies of study counts
for version in Releases:
    CNVGreaterTargetscan.append(CompTargetscan[version][0] / sum(CompTargetscan[version]))
    CNVLowerTargetscan.append(CompTargetscan[version][1] / sum(CompTargetscan[version]))
    CNVNoDiffTargetscan.append(CompTargetscan[version][2] / sum(CompTargetscan[version]))
    CNVGreaterMiranda.append(CompMiranda[version][0] / sum(CompMiranda[version]))
    CNVLowerMiranda.append(CompMiranda[version][1] / sum(CompMiranda[version]))
    CNVNoDiffMiranda.append(CompMiranda[version][2] / sum(CompMiranda[version]))

       
# make list of data
TargetscanData = [CNVGreaterTargetscan, CNVLowerTargetscan, CNVNoDiffTargetscan]
MirandaData = [CNVGreaterMiranda, CNVLowerMiranda, CNVNoDiffMiranda]       



# create figure
fig = plt.figure(1, figsize = (6, 3))

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
    # make a list with added values for greater and lower
    added = []
    for i in range(len(greater)):
        added.append(greater[i] + lower[i])
    # Create a bar plot for proportions of studies with CNV lower on top of CNV greater
    ax.bar([0, 0.4, 0.8, 1.2], greater, width = 0.3, label = 'CNV > non-CNV', color= '#ef8a62')
    # Create a bar plot for proportions of studies with CNV lower on top of CNV greater
    ax.bar([0, 0.4, 0.8, 1.2], lower, width = 0.3, bottom = greater, label = 'CNV < non-CNV', color = '#67a9cf')
    # Create a bar plot with proportions of studies with no difference between CNV and non-CNV genes on top of the other 2 bars
    ax.bar([0, 0.4, 0.8, 1.2], nodiff, width = 0.3, bottom = added, label = 'No difference', color= '#f7f7f7')
    
    # write title
    ax.set_title(Title + '\n', size = 10)
    # set font for all text in figure
    FigFont = {'fontname':'Arial'}   
    # write label for y axis
    ax.set_ylabel('Frequency of studies in DGV', color = 'black', size = 10, ha = 'center', **FigFont)
    # write label for x axis
    ax.set_ylabel('DGV release', color = 'black', size = 10, ha = 'center', **FigFont)
    # add X axis tick labels
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
ax1 = CreateAx(2, 1, 1, TargetscanData, fig, 'TargetScan', labelnames, [0.15, 0.55, 0.95, 1.35])
ax2 = CreateAx(2, 1, 2, MirandaData, fig, 'miRanda', labelnames, [0.15, 0.55, 0.95, 1.35])

# add subplot label
ax1.text(-0.5, 1.1, 'A', horizontalalignment = 'center',
         verticalalignment = 'center', color = 'black', size = 10)
ax1.text(2, 1.1, 'B', horizontalalignment = 'center',
         verticalalignment = 'center', color = 'black', size = 10)

# add legend
N = mpatches.Patch(facecolor = '#f7f7f7' , edgecolor = 'black', linewidth = 1, label= 'No diff.')
G = mpatches.Patch(facecolor = '#ef8a62' , edgecolor = 'black', linewidth = 1, label= 'CNV greater')
L = mpatches.Patch(facecolor = '#67a9cf' , edgecolor = 'black', linewidth = 1, label= 'CNV lower')
ax1.legend(handles = [G, L, N], loc = (0, 1), fontsize = 8, frameon = False, ncol = 3)

# make sure subplots do not overlap
plt.tight_layout()

## get outputfile
#outputfile = 'truc' + '_' + domain + '_' + chromos + '_' + cnv_length 
#print(outputfile)

# save figure
fig.savefig('truc.pdf', bbox_inches = 'tight')
       
