# -*- coding: utf-8 -*-
"""
Created on Fri Jul 22 17:40:04 2016

@author: RJovelin
"""

# use this script to plot the number of short 3'UTR genes in CNV and non-CNVs for release of the DGV

# usage PlotNumberShort3UTRGenes.py [options]
# [7/15]: 3'UTR length, genes with 3'UTR length < 7bp or < 15bp are considered short 3'UTR genes

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
# alternatives
# chromos = 'all_chromos'
# cnv_length = 'CNV_greater_1Kb'

# get the minimum 3'UTR length
L = int(sys.argv[1])
assert L in [15, 7], 'minimum 3UTR length is not correct'
print(L)

# Count the number of short 3'UTR genes in CNV and non-CNV in each release
# create a dict {release: [N short 3'UTR, N short 3UTR CNV genes, N short 3UTR non-CNV genes]}
ShortGenes = {}

# get UTR file
UTR_file = 'H_sapiens_3UTR_length_' + chromos + '.txt'
print(UTR_file)

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
    # sort genes based on 3' UTR length
    UTR_length = sort_genes_3UTR_length(UTR_file, L)
    print(len(UTR_length))
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

# create figure
fig = plt.figure(1, figsize = (4.3,2.56))

# add axe to fig
ax = fig.add_subplot(1, 1, 1)

# Set the bar width
bar_width = 0.5

# set positions of the x-axis ticks
xtickpos = [0, 0.7, 1.4, 2.1]

# Create a bar plot for cnv genes
ax.bar(xtickpos, cnv_genes, width=bar_width, label = 'CNV', color= '#ef8a62')
# Create a bar plot for non_cnv genes on top of cnv_genes
ax.bar(xtickpos, non_cnv_genes, width=bar_width, bottom= cnv_genes, label = 'non-CNV', color = '#67a9cf')

# set font for all text in figure
FigFont = {'fontname':'Arial'} 

# add labels to x axis ticks
plt.xticks(xtickpos, labelnames, **FigFont)

# set axis labels
ax.set_ylabel('Number of genes\nwith short 3\'UTR', size = 10, ha = 'center', color = 'black', **FigFont)

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

# get outputfile
if L == 7:
    outputfile = 'PlotShort3UTRCountsDGV_7bp_' + cnv_length + '_' + chromos
elif L == 15:
    outputfile = 'PlotShort3UTRCountsDGV_15bp_' + cnv_length + '_' + chromos
print(outputfile)
  
# save figure
fig.savefig(outputfile + '.eps', bbox_inches = 'tight')
