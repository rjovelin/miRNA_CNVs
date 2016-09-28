# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 14:10:59 2016

@author: RJovelin
"""


# use this script to compare the number target sites per nucleotide between CNV and non-CNV genes
# for single human populations 

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
import copy
# import custom modules
from CNV_miRNAs import *


# usage PlotTargetsSingleHumanPop.py [options]
# [3UTR/5UTR/CDS] choose the region to analyse
# [pdf/ai/png]: save as eps if no argument is provided or as the format indicated by the argument

# get the region to consider to predict target sites [3UTR or 5UTr or CDS]
domain = sys.argv[1]
print(domain)
# check the format of the figure (eps by default or pdf if specified)
if len(sys.argv) == 2:
    extension = '.eps'
elif len(sys.argv) == 3:
    extension = sys.argv[2]
    assert extension in ['pdf', 'ai', 'png']
    extension = '.' + extension

# keep genes on assembled chromsomes
keep_valid_chromos = True
chromos = 'valid_chromos'
# use all CNVs
cnv_length = 'CNV_all_length'

# get the CNV files with single human population references
CNV_file = 'GRCh37_hg19_variants_2016-05-15.txt'

# get the file with 3'UTR length
UTR_file = 'H_sapiens_3UTR_length_' + chromos + '.txt'
print(UTR_file)
# get targetscan sequence input file
targetscan_seq_input_file = 'H_sapiens_' + domain + '_' + chromos + '_targetscan.txt'
print(targetscan_seq_input_file)
# get the outputfiles with predicted target sites for targetscan and miranda
TargetscanFile = 'H_sapiens_' + domain + '_' + chromos + '_predicted_sites_targetscan.txt'    
MirandaFile = 'H_sapiens_' + domain + '_' + chromos + '_predicted_sites_miranda.txt'
print(TargetscanFile, MirandaFile)
# make a dictionary with {gene :[targets, seq_length, normalized_targets]}
TargetscanTargets = parse_targetscan_output(targetscan_seq_input_file, TargetscanFile , 'all')
MirandaTargets = parse_miranda_output(targetscan_seq_input_file, MirandaFile, 'all')
print('# targetscan sites', len(TargetscanTargets))
print('# miranda sites', len(MirandaTargets))


# make a dictionary {reference: pubmedid} for studies reported in the CGV
references = get_DGV_references(CNV_file)

# make a list of pubmed IDs for single population studies
pubmed = [['Suktitipat_et_al_2014', '25118596'], ['John_et_al_2014', '26484159'], ['Thareja_et_al_2015', '25765185'], ['Alsmadi_et_al_2014', '24896259']]

# check that reference names correspond to the pubmed IDs
infile = open(CNV_file)
for line in infile:
    for i in range(len(pubmed)):
        if pubmed[i][0] in line:
            assert pubmed[i][1] in line, 'study ref does not correspond to expected pubmed ID'
infile.close()


# create a set of CNV genes for each study {study: {gene1, gene2...}}
StudiesCNV = {}
for i in range(len(pubmed)):
    # get the set of CNV genes for that study
    CNV_genes = get_human_CNV_genes_single_study(CNV_file, pubmed[i][0], 'all')
    print(pubmed[i][0], len(CNV_genes))
    # populate dict
    StudiesCNV[pubmed[i][0]] = CNV_genes
print('extracted CNV genes for each study')

# get synonym names for all genes {gene name : [list of synonyms]}
synonyms = get_synonyms('H_sapiens.gff3')
print('got gene synonyms', len(synonyms))    
# get the CDS sequences of the longest mRNAs for each gene {gene : sequence}
CDS_seq = extract_CDS_sequences('H_sapiens.gff3', 'H_sapiens_genome.txt', 'H_sapiens_valid_chromos.txt', keep_valid_chromos)
print('extracted CDS sequences', len(CDS_seq))  


# create a dict {study: {gene: CNV status}}    
CNV_status = {}
for study in StudiesCNV:
    # initialize dict
    CNV_status[study] = {}
    # loop over genes in CDS seq
    for gene in CDS_seq:
        # set boolean
        is_cnv = False
        # ask if gene in CNV genes
        if gene in StudiesCNV[study] or gene.upper() in StudiesCNV[study]:
            # gene is CNV, add gene and status to inner dict
            CNV_status[study][gene] = 'CNV'
        else:
            # ask if any of the gene synonyms are in CNV genes
            for name in synonyms[gene]:
                # check if gene in CNV
                if name in StudiesCNV[study] or name.upper() in StudiesCNV[study]:
                    # update boolean variable
                    is_cnv = True
            # check if gene in CNV
            if is_cnv == True:
                CNV_status[study][gene] = 'CNV'
            elif is_cnv == False:
                CNV_status[study][gene] = 'not_CNV'
print('sorted genes according to CNV status')    
    
    
# use this function to add the cnv status to each gene with predicted targets
def CombineTargetsCNV(targets, CNV_status):
    '''
    (dict, dict) - > dict
    Take a dictionary with predicted targets for each gene and a dictionary with
    CNV status for each gene in each study and a return a dictionary with targets
    and CNV status for each gene in each study     
    '''
    CNVTargets = {}
    for study in CNV_status:
        # initialize inner dict
        CNVTargets[study] = {}
        # loop over genes with CNV status for given study
        for gene in CNV_status[study]:
            # populate with list of targets
            if gene in targets:
                CNVTargets[study][gene] = copy.deepcopy(targets[gene])
                # add CNV status
                CNVTargets[study][gene].append(CNV_status[study][gene])
    return CNVTargets


# use this function count the number of CNV genes and non-CNV genes in each study
def CountGenes(CNVTargets):
    '''
    (dict) -> dict
    Take a dictionary with CNV status and predicted targets for each study
    and return a dict with the counts of CNV and non_CNV genes in each study
    '''
    CNVNum = {}
    for study in CNVTargets:
        CNVNum[study] = [0, 0]
        for gene in CNVTargets[study]:
            if CNVTargets[study][gene][-1] == 'CNV':
                CNVNum[study][0] += 1
            elif CNVTargets[study][gene][-1] == 'not_CNV':
                CNVNum[study][1] += 1
    return CNVNum


# use this function to combine the number of predicted targets in arrays
def MakeTargetArray(CNVTargets):
    '''
    (dict) -> dict
    Take a dict with the predicted targets and CNV status for each gene in each study
    and return a dict with lists of targets for CNV and non-non-CNV genes in each study    
    '''
    # create a dict of {study : [[CNV], [non-CNV]]}
    CNVData = {} 
    # loop over studies, get the number of normalized sites for CNV and non-CNV genes
    for study in CNVTargets:
        # initialise list value
        CNVData[study] = [[], []]
        # populate inner lists with number of miRNA target sites per nucleotide
        for gene in CNVTargets[study]:
            if CNVTargets[study][gene][-1] == 'CNV':
                CNVData[study][0].append(CNVTargets[study][gene][2])
            elif CNVTargets[study][gene][-1] == 'not_CNV':
                CNVData[study][1].append(CNVTargets[study][gene][2])
    return CNVData
    

    

# make a dictionary {study: {gene: [targets, seq_length, normalized_targets, CNV_status]}}
TargetscanCNVTargets = CombineTargetsCNV(TargetscanTargets, CNV_status)
MirandaCNVTargets = CombineTargetsCNV(MirandaCNVTargets, CNV_status)
print('combined targets and CNV status')

# count CNV and non-CNV genes in each study
CNVNumTargetscan = CountGenes(TargetscanCNVTargets)
CNVNumMiranda = CountGenes(MirandaCNVTargets)
print('got CNV gene counts for each study')
for study in CNVNumTargetscan:
    print(study, CNVNumTargetscan[study][0], CNVNumTargetscan[study][1])
    print(study, CNVNumMiranda[study][0], CNVNumMiranda[study][1])

# create a dict of {study : [[CNV], [non-CNV]]}
CNVDataTargetscan = MakeTargetArray(TargetscanCNVTargets)
CNVDataMiranda = MakeTargetArray(MirandaCNVTargets) 
print('generated lists of target sites for CNV and non-CNV genes')






# make a list of all data for each study
AllDataTargetscan, AllDataMiranda = [], []
# make a list of species names to loop from
StudyNames = ['Suktitipat_et_al_2014', 'Alsmadi_et_al_2014', 'John_et_al_2014', 'Thareja_et_al_2015']
Populations = ['$Thai^a$', '$Kuwaiti^b$', '$Kuwaiti^c$', '$Kuwaiti^d$']

# loop over study in studies list and populate data lists, keeping the same order for targetscan and miranda
for study in StudyNames:
    # append list of target sites for CNV genes
    AllDataTargetscan.append(CNVDataTargetscan[study][0])
    AllDataMiranda.append(CNVDataMiranda[study][0])    
    # append list of target sites for non-CNV genes
    AllDataTargetscan.append(CNVDataTargetscan[study][1])
    AllDataMiranda.append(CNVDataMiranda[study][1])
print('data consolidated in array')


# boostrap CNV and non-CNV genes to compare miRNA targets
# create a dictionary {study: {CNV_status: {num: gene}}} 
ToSampleFromTargetscan, ToSampleFromMiranda = {}, {}
# loop over studies
for study in TargetscanCNVTargets:
    # initialize dict
    ToSampleFromTargetscan[study] = {}
    ToSampleFromTargetscan[study]['CNV'], ToSampleFromTargetscan[study]['not_CNV'] = {}, {}
    # initialize counters for CNV and non-CNV genes
    i, j = 0, 0
    # loop over genes for that study
    for gene in TargetscanCNVTargets[study]:
        if TargetscanCNVTargets[study][gene][-1] == 'CNV':
            # populate dict with num : gene pair and update counter
            ToSampleFromTargetscan[study]['CNV'][i] = gene
            i += 1
        elif TargetscanCNVTargets[study][gene][-1] == 'not_CNV':
            # populate dict with num : gene pair and update counter
            ToSampleFromTargetscan[study]['not_CNV'][j] = gene
            j += 1
for study in MirandaCNVTargets:
    # initialize dict
    ToSampleFromMiranda[study] = {}
    ToSampleFromMiranda[study]['CNV'], ToSampleFromMiranda[study]['not_CNV'] = {}, {}
    # initialize counters for CNV and non-CNV genes
    i, j = 0, 0
    # loop over genes for that study
    for gene in MirandaCNVTargets[study]:
        if MirandaCNVTargets[study][gene][-1] == 'CNV':
            # populate dict with num : gene pair and update counter
            ToSampleFromMiranda[study]['CNV'][i] = gene
            i += 1
        elif MirandaCNVTargets[study][gene][-1] == 'not_CNV':
            # populate dict with num : gene pair and update counter
            ToSampleFromMiranda[study]['not_CNV'][j] = gene
            j += 1
print('assigned numbers to gene for sampling')

# check that all genes have been assigned to a number
for study in ToSampleFromTargetscan:
    assert len(ToSampleFromTargetscan[study]['CNV']) == CNVNumTargetscan[study][0], 'targetscan CNV genes counts do not match'
    assert len(ToSampleFromTargetscan[study]['not_CNV']) == CNVNumTargetscan[study][1], 'targetscan non-CNV genes counts do not match'
    assert len(ToSampleFromMiranda[study]['CNV']) == CNVNumMiranda[study][0], 'miranda CNV genes counts do not match'
    assert len(ToSampleFromMiranda[study]['not_CNV']) == CNVNumMiranda[study][1], 'miranda non-CNV genes counts do not match'















# create a dict for each study with a list with numbers of each different outcomes when comparing targets in CNV and non-CNV genes
# {study: [# replicates CNV > non-CNV, # replicates CNV < non-CNV, # replicates no differences]}
BootStrap = {}
# initialize list values
for study in ToSampleFrom:
    BootStrap[study] = [0, 0, 0]

# loop over studies in dict to sample from
for study in ToSampleFrom:
    print('bootstraping', study)
    # set number of replicates
    replicates = 10000
    while replicates != 0:
        # make list of targets for CNV and non-CNV genes
        repCNVtargets, repNonCNVtargets = [], []
        # draw 800 CNV genes and 800 non-CNV genes with replacement
        for i in range(800):
            # draw a random CNV gene
            j = random.randint(0, len(ToSampleFrom[study]['CNV']) - 1)
            k = random.randint(0, len(ToSampleFrom[study]['not_CNV']) - 1)
            # get the corresponding genes
            gene1 = ToSampleFrom[study]['CNV'][j]
            gene2 = ToSampleFrom[study]['not_CNV'][k]            
            # get the the number of targets for these 2 genes
            assert CNVTargets[study][gene1][-1] == 'CNV', 'random gene should be CNV'
            assert CNVTargets[study][gene2][-1] == 'not_CNV', 'random gene should be non-CNV'
            repCNVtargets.append(CNVTargets[study][gene1][2])
            repNonCNVtargets.append(CNVTargets[study][gene2][2])
        # make sure that the correct numbers of genes is drawn
        assert len(repCNVtargets) == 800, '800 CNV genes should be drawn'
        assert len(repNonCNVtargets) == 800, '800 non-CNV genes should be drawn'
        # compare CNV and non-CNV genes
        Pval = stats.ranksums(repCNVtargets, repNonCNVtargets)[1]
        # check significance
        if Pval >= 0.05:
            # difference is not significance
            BootStrap[study][2] += 1
        elif Pval < 0.05:
            # difference is significance, check if CNV genes have a greater number of targets
            if np.mean(repCNVtargets) > np.mean(repNonCNVtargets):
                BootStrap[study][0] += 1
            elif np.mean(repCNVtargets) < np.mean(repNonCNVtargets):
                BootStrap[study][1] += 1
        # update replicate number
        replicates -= 1
print('done with boostraping')        

      
# create parallel list of proportions 
Greater, Lower, Nodiff = [], [] ,[]
for study in StudyNames:
    Greater.append(BootStrap[study][0] / sum(BootStrap[study]))
    Lower.append(BootStrap[study][1] / sum(BootStrap[study]))
    Nodiff.append(BootStrap[study][2] / sum(BootStrap[study]))

# create a list with all the proportion lists
Proportions = [Greater, Lower, Nodiff]

# print results to screen
for study in BootStrap:
    print(BootStrap[study])



# create figure
fig = plt.figure(1, figsize = (8, 3))

# create a function to format the subplots
def CreateAx(Columns, Rows, Position, Data, figure, Title, YMax, LabelNames, XScale, GraphType):
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
    
    if GraphType == 'box':
        # create a list of positions for the box plot    
        BoxPositions = [0, 0.4, 0.9, 1.3, 1.8, 2.2, 2.7, 3.1]
        # use a boxplot
        bp = ax.boxplot(Data, showmeans = True, showfliers = False, widths = 0.3,
                        positions = BoxPositions, patch_artist = True) 
    
        # color CNV and non-CNV boxes differently
        i = 0    
        # change box, whisker color to black
        for box in bp['boxes']:
            # change line color
            box.set(color = 'black')
            if i % 2 == 0:
                # CNV data, color box in grey
                box.set(facecolor = '#a6cee3')
            else:
                box.set(facecolor = '#b2df8a')
            i += 1
        # change whisker color to black
        for wk in bp['whiskers']:
            wk.set(color = 'black', linestyle = '-')
        # change color of the caps
        for cap in bp['caps']:
            cap.set(color = 'black')
        # change the color and line width of the medians
        for median in bp['medians']:
            median.set(color = 'black')
        # change the mean marker and marker
        for mean in bp['means']:
            mean.set(marker = 'o', markeredgecolor = 'black', markerfacecolor = 'black', markersize = 3)
    
        # add legend
        C = mpatches.Patch(facecolor = '#a6cee3' , edgecolor = 'black', linewidth = 1, label= 'CNV')
        N = mpatches.Patch(facecolor = '#b2df8a' , edgecolor = 'black', linewidth = 1, label= 'non-CNV')
        plt.legend(handles = [C, N], loc = (0, 1), fontsize = 8, frameon = False, ncol = 2)
    
  
    elif GraphType == 'bar':
        # get the list of proportions 
        greater, lower, nodiff = Data[0], Data[1], Data[2]
        # make a list with added values for nodiff and greater
        added = []
        for i in range(len(greater)):
            added.append(nodiff[i] + greater[i])
        # Create a bar plot for proportions of replicates with CNV no diff on top of CNV lower
        ax.bar([0, 0.4, 0.8, 1.2], nodiff, width = 0.3, label = 'No difference', color= '#f7f7f7')
        # Create a bar plot for proportions of replicates with CNV greater on top of no diff
        ax.bar([0, 0.4, 0.8, 1.2], greater, width = 0.3, bottom = nodiff, label = 'CNV > non-CNV', color= '#ef8a62')
        # Create a bar plot for proportions of replicates with CNV lower on top of CNV greater
        ax.bar([0, 0.4, 0.8, 1.2], lower, width = 0.3, bottom= added, label = 'CNV < non-CNV', color = '#67a9cf')
 
        # add legend
        N = mpatches.Patch(facecolor = '#f7f7f7' , edgecolor = 'black', linewidth = 1, label= 'No diff.')
        G = mpatches.Patch(facecolor = '#ef8a62' , edgecolor = 'black', linewidth = 1, label= 'CNV greater')
        L = mpatches.Patch(facecolor = '#67a9cf' , edgecolor = 'black', linewidth = 1, label= 'CNV lower')
        plt.legend(handles = [N, G, L], loc = (0, 1), fontsize = 8, frameon = False, ncol = 3)

    # write title   
    ax.set_title(Title + '\n\n', size = 8)
    
    # set font for all text in figure
    FigFont = {'fontname':'Arial'}   
    
    # write label for y axis
    if GraphType == 'box':
        ax.set_ylabel('Normalized number of miRNA\nsites per gene', color = 'black',  size = 8, ha = 'center', **FigFont)
    elif GraphType == 'bar':
        ax.set_ylabel('Proportion of replicates', color = 'black', size = 8, ha = 'center', **FigFont)

    # write label for x axis
    plt.xticks(XScale, LabelNames, ha = 'center', fontsize = 8, **FigFont)

    # add a range for the Y amd X axes
    if GraphType == 'box':
        plt.ylim([0, YMax])
        plt.xlim([-0.25, 3.35])
    
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
        labelsize = 8,
        direction = 'out') # ticks are outside the frame when bottom = 'on'  
      
    # Set the tick labels font name
    for label in ax.get_yticklabels():
        label.set_fontname('Arial')
    
    # create a margin around the x axis
    plt.margins(0.05)
    
    return ax      



# plot boxplots for predictor in 1st subplot
# get title based on predictor algorithm
if predictor == 'targetscan':
    figtitle = 'TargetScan'
elif predictor == 'miranda':
    figtitle = 'miRanda'
ax1 = CreateAx(2, 1, 1, AllData, fig, figtitle, 0.45, Populations, [0.2, 1.1, 2, 2.9], 'box')
# plot bars in 2nd subplot
ax2 = CreateAx(2, 1, 2, Proportions, fig, figtitle, 0.45, Populations, [0.15, 0.55, 0.95, 1.35], 'bar')


# annotate Graph with significance level
Pvalues = []
for study in StudyNames:
    # compare mirna targets between CNV and non-CNV genes
    P = stats.ranksums(CNVData[study][0], CNVData[study][1])[1]
    if P >= 0.05:
        Pvalues.append('')
    elif P < 0.05 and P >= 0.01:
        Pvalues.append('*')
    elif P < 0.01 and P >= 0.001:
        Pvalues.append('**')
    elif P < 0.001:
        Pvalues.append('***')
# create list of Y and X positions to annotate figure with significance level
Ypos = [0.42, 0.42, 0.42, 0.42]
Xpos = [0.2, 1.1, 2, 2.9]
for i in range(len(Pvalues)):
    ax1.text(Xpos[i], Ypos[i], Pvalues[i], horizontalalignment = 'center',
             verticalalignment = 'center', color = 'black', size = 8)


# add subplot label
ax1.text(-1, 0.48, 'A', horizontalalignment = 'center',
         verticalalignment = 'center', color = 'black', size = 10)
ax1.text(3.5, 0.48, 'B', horizontalalignment = 'center',
         verticalalignment = 'center', color = 'black', size = 10)


# make sure subplots do not overlap
plt.tight_layout()

# get outputfile
outputfile = 'PlotSinglePops_' + predictor + '_' + domain + '_' + chromos + '_' + cnv_length 
print(outputfile)

# save figure
fig.savefig(outputfile + extension, bbox_inches = 'tight')
       
