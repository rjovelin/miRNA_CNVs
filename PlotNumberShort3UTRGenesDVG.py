# -*- coding: utf-8 -*-
"""
Created on Fri Jul 22 17:40:04 2016

@author: RJovelin
"""

# use this script to plot the number of short 3'UTR genes in CNV and non-CNVs for release of the DGV
# and the ratio of CNV / non-CNV genes per study for each version of the DGV

# usage PlotNumberShort3UTRGenes.py [options]
# [7/15]: 3'UTR length, genes with 3'UTR length < 7bp or < 15bp are considered short 3'UTR genes
# [/pdf]: save as eps if no argument is provided or as PDF if pdf is given in the command

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
# import custom modules
from CNV_miRNAs import *


# use all chromos (including unplaced, unlocated, and MT) or only valid chromos 
# note: only CNVs on valid chromos are reported in DGV, so if all chromos are
# used it may introduce a bias by calling non CNV genes genes that cannot be detected

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

# get the minimum 3'UTR length
L = int(sys.argv[1])
assert L in [15, 7], 'minimum 3UTR length is not correct'
print(L)
# check the format of the figure (eps by default or pdf if specified)
if len(sys.argv) == 2:
    extension = '.eps'
elif len(sys.argv) == 3:
    assert sys.argv[2] == 'pdf'
    extension = '.pdf'


# Count the number of short 3'UTR genes in CNV and non-CNV in each release
# create a dict {release: [N short 3'UTR, N short 3UTR CNV genes, N short 3UTR non-CNV genes]}
ShortGenes = {}

# get UTR file
UTR_file = 'H_sapiens_3UTR_length_' + chromos + '.txt'
print(UTR_file)

# sort genes based on 3' UTR length
UTR_length = sort_genes_3UTR_length(UTR_file, L)
print(len(UTR_length))

# make a list of CNV files
CNV_files = ['H_sapiens_GRCh37_2013-05_' + cnv_length + '_' + chromos + '.txt',
             'H_sapiens_GRCh37_2013-07_' + cnv_length + '_' + chromos + '.txt',
             'H_sapiens_GRCh37_2014_' + cnv_length + '_' + chromos + '.txt',
             'H_sapiens_GRCh37_2015_' + cnv_length + '_' + chromos + '.txt']

# loop over CNV files
for filename in CNV_files:
    print(filename)
    # sort genes based on CNV status
    CNV_status = sort_genes_CNV_status(filename)
    print(len(CNV_status))
    # get release version
    release_version = filename[filename.index('GRCh37'): filename.index('_CNV')]
    print(release_version)
    # count total number of genes with short 3'UTR
    total_short = 0
    # count CNV genes with short UTR
    cnv_short = 0
    # count non-CNV genes with short UTR
    non_cnv_short = 0    
    # loop over genes in UTR_length
    for gene in UTR_length:
        if UTR_length[gene] == 'short':
            total_short += 1
            # check CNV status
            if CNV_status[gene] == 'CNV':
                cnv_short += 1
            elif CNV_status[gene] == 'not_CNV':
                non_cnv_short += 1
            
    # check that numbers add up
    assert total_short == cnv_short + non_cnv_short, 'sum cnv and non-cnv short is not equal to total short'
    # populare dict
    ShortGenes[release_version] = [total_short, cnv_short, non_cnv_short]
print('got short gene counts for each species')


# plot the number of genes with short 3'UTR in each release of the DVG as a bar graph

Releases = ['GRCh37_2013-05', 'GRCh37_2013-07', 'GRCh37_2014', 'GRCh37_2015']
labelnames = ['2013a', '2013b', '2014', '2015']
    
# create list of counts for cnv genes and non-CNV genes parallel to labelnames list
cnv_genes, non_cnv_genes = [], []
for i in Releases:
    cnv_genes.append(ShortGenes[i][1])
    non_cnv_genes.append(ShortGenes[i][2])

# make a list of cnv data
CNVData = [cnv_genes, non_cnv_genes]


# count short 3'UTR genes for each study of each release of the DGV

# make a list of DGV files
DGVFiles = ['GRCh37_hg19_variants_2013-05-31.txt', 'GRCh37_hg19_variants_2013-07-23.txt',
            'GRCh37_hg19_variants_2014-10-16.txt', 'GRCh37_hg19_variants_2015-07-23.txt']

# create a dict {release: {study: [pubmedID, total, CNV, non-CNV]}}
StudiesShortGenes = {}

# create a dict {release: {study: {set of cnv genes}}}
StudiesCNVGenes = {}

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

# get synonym names for all genes {gene name : [list of synonyms]}
synonyms = get_synonyms('H_sapiens.gff3')
print('got synonymous names', len(synonyms))

# get the CDS sequences of the longest mRNAs for each gene {gene : sequence}
CDS_seq = extract_CDS_sequences('H_sapiens.gff3', 'H_sapiens_genome.txt', 'H_sapiens_valid_chromos.txt', keep_valid_chromos)
print('extracted CDS sequences', len(CDS_seq))  


# get the CNV genes for each study of each release of the DGV
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
    
# create a dict with ratio CNV / non_CNV genes for each study and release
# {release: [list of ratio CNV / non-CNV genes]}
StudiesRatio = {}
for release in CNV_status:
    StudiesRatio[release] = []
    for study in CNV_status[release]:
        cnvcount, noncnvcount = 0, 0
        for gene in CNV_status[release][study]:
            if CNV_status[release][study][gene] == 'CNV':
                cnvcount += 1
            elif CNV_status[release][study][gene] == 'not_CNV':
                noncnvcount += 1
        # populate dict
        StudiesRatio[release].append(cnvcount / noncnvcount)
print('got the CNV / non-CNV gene ratio for each study of each release')

# get the maximum ratio value
ratio = []
for i in StudiesRatio:
    for j in StudiesRatio[i]:
        ratio.append(j)
MaxRatio = max(ratio)
print(MaxRatio, (((MaxRatio * 100) // 10) + 1))


# plot the number of studies with ratio of CNV / non-CNV short 3'UTR genes
# for each release of the DGV

#create a dict with version and list of counts {release: [n1, n2, etc]}
Ratio = {}
for release in StudiesRatio:
    # create a list of 0
    counts = [0] * int(((MaxRatio * 100) // 10) + 1)
    # get the index in list where value should go
    for freq in StudiesRatio[release]:
        pos = int((freq * 100) // 10)
        counts[pos] += 1
    # populate dict
    Ratio[release] = counts

    
# make a list of frequency data
RatioData = [Ratio, Releases, labelnames]

# create figure
fig = plt.figure(1, figsize = (7, 2.5))


# create a function to format the subplots
def CreateAx(Columns, Rows, Position, Data, figure, YMax, LabelNames, XScale, GraphType):
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
    
    
    if GraphType == 'bar':
        # set positions of the x-axis ticks
        xtickpos = [0, 0.7, 1.4, 2.1]
        # parse list Data
        cnv_genes, non_cnv_genes = Data[0], Data[1]
        # Create a bar plot for cnv genes
        ax.bar(xtickpos, cnv_genes, width= 0.5, label = 'CNV', color= '#ef8a62')
        # Create a bar plot for non_cnv genes on top of cnv_genes
        ax.bar(xtickpos, non_cnv_genes, width= 0.5, bottom= cnv_genes, label = 'non-CNV', color = '#67a9cf')
    
        # add legend
        C = mpatches.Patch(facecolor = '#ef8a62', edgecolor = 'black', linewidth = 1, label= 'CNV')
        N = mpatches.Patch(facecolor = '#67a9cf', edgecolor = 'black', linewidth = 1, label= 'non-CNV')
        plt.legend(handles = [C, N], loc = (0, 1), fontsize = 8, frameon = False, ncol = 2)
    
    elif GraphType == 'line':
        # create a dict to build the legend
        Graphics = {}
        # parse list data
        Ratio, Releases, labelnames = Data[0], Data[1], Data[2]
        for i in range(len(Releases)):
            if '2013-05' in Releases[i]:
                graph = ax.plot([j + 0.5 for j in range(len(Ratio[Releases[i]]))], Ratio[Releases[i]], linestyle = '-', color = '#b3cde3', marker = 'o', markersize = 4, markeredgewidth = 1, markerfacecolor = '#b3cde3', markeredgecolor = '#b3cde3', lw = 1.5, label = LabelNames[i]) 
            elif '2013-07' in Releases[i]:
                graph = ax.plot([j + 0.5 for j in range(len(Ratio[Releases[i]]))], Ratio[Releases[i]], linestyle = '-', color = '#8c96c6', marker = 'o', markersize = 4, markeredgewidth = 1, markerfacecolor = '#8c96c6', markeredgecolor = '#8c96c6', lw = 1.5, label = LabelNames[i]) 
            elif '2014' in Releases[i]:
                graph = ax.plot([j + 0.5 for j in range(len(Ratio[Releases[i]]))], Ratio[Releases[i]], linestyle = '-', color = '#8856a7', marker = 'o', markersize = 4, markeredgewidth = 1, markerfacecolor = '#8856a7', markeredgecolor = '#8856a7', lw = 1.5, label = LabelNames[i]) 
            elif '2015' in Releases[i]:
                graph = ax.plot([j + 0.5 for j in range(len(Ratio[Releases[i]]))], Ratio[Releases[i]], linestyle = '-', color = '#810f7c', marker = 'o', markersize = 4, markeredgewidth = 1, markerfacecolor = '#810f7c', markeredgecolor = '#810f7c', lw = 1.5, label = LabelNames[i]) 
            # populate dict
            Graphics[labelnames[i]] = graph

        # add legend
        # add lines
        lns = Graphics[labelnames[0]]
        for i in labelnames[1:]:
            lns += Graphics[i]
        # get labels
        labs = [i for i in labelnames]
        ax.legend(lns, labs, loc = 1, fontsize = 8, frameon = False, numpoints = 1)    

    # set font for all text in figure
    FigFont = {'fontname':'Arial'}   
    
    # write label for x and y axis
    if GraphType == 'bar':
        # set axis labels
        ax.set_ylabel('Number of genes\nwith short 3\'UTR', size = 10, ha = 'center', color = 'black', **FigFont)
        ax.set_xlabel('DGV releases', fontsize = 10, **FigFont)
    elif GraphType == 'line':
        ax.set_ylabel('Number of studies in DGV', size = 10, ha = 'center', color = 'black', **FigFont)
        ax.set_xlabel('CNV genes / non-CNV genes ratio', fontsize = 10, **FigFont)

    # write label for x axis
    plt.xticks(XScale, LabelNames, ha = 'center', fontsize = 8, **FigFont)
    
    # limit y axis range
    if GraphType == 'bar':
        plt.ylim([0, YMax])
    elif GraphType == 'line':
        plt.ylim([-1, YMax])
        
    
    # do not show lines around figure  
    ax.spines["top"].set_visible(False)    
    ax.spines["bottom"].set_visible(True)    
    ax.spines["right"].set_visible(False)    
    ax.spines["left"].set_visible(False)  
    # offset the spines
    for spine in ax.spines.values():
        spine.set_position(('outward', 5))

    # add a light grey horizontal grid to the plot, semi-transparent, 
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)  
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


# plot the number of short 3' UTR genes
ax1 = CreateAx(2, 1, 1, CNVData, fig, 700, labelnames, [0.25, 0.95, 1.65, 2.35], 'bar')
# plot the ratio of CNV / non-CNV genes
subplot2xticks = []
for i in range(len(Ratio['GRCh37_2013-05']) + 1):
    if i % 2 == 0:
        subplot2xticks.append(str(i/10))
    else:
        subplot2xticks.append('')
ax2 = CreateAx(2, 1, 2, RatioData, fig, 50, subplot2xticks, [i for i in range(len(Ratio['GRCh37_2013-05']) + 1)], 'line')

# make sure subplots do not overlap
plt.tight_layout()  
   
   
# add subplot label
ax1.text(-1, 800, 'A', horizontalalignment = 'center',
         verticalalignment = 'center', color = 'black', size = 10)
ax1.text(3, 800, 'B', horizontalalignment = 'center',
         verticalalignment = 'center', color = 'black', size = 10)   
   
# get outputfile
if L == 7:
    outputfile = 'PlotShort3UTRCountsDGV_7bp_' + cnv_length + '_' + chromos
elif L == 15:
    outputfile = 'PlotShort3UTRCountsDGV_15bp_' + cnv_length + '_' + chromos
print(outputfile)

# save figure
fig.savefig(outputfile + extension, bbox_inches = 'tight')
