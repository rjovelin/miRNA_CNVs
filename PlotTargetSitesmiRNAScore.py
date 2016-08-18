# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 15:55:09 2016

@author: RJovelin
"""



# use this script to compare target sites between human CNV and non-CNv genes
# assigning different scores to miRNAs according to their expression level

# usage PlotTargetSitesmiRNAScore.py [options]
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


fastafile = 'H_sapiens_miRBaseMatureAccession.txt'
# create a dict mature name : score pairs
scores = TargetScore(fastafile, 1000, 'Homo_sapiens', miRBaseFile = 'miRNA.dat', ExpressionFile = 'mirna_read_count.txt')
print('match scores to mature miRNAs', len(scores))

# get the seq input file
seq_input_file = 'H_sapiens_' + domain + '_' + chromos + '_targetscan.txt'
# get the predicted targets output file
predicted_targets = 'H_sapiens_' + domain + '_' + chromos + '_predicted_sites_miranda.txt'
# use DGV 2015 release 
CNV_file = 'H_sapiens_GRCh37_2015_CNV_all_length_valid_chromos.txt'

# record the number of miranda target sites for each gene weighted by the mirna expression score
#  {gene: [N_targets, Sequence_length, N_targets_normalized, CNV_status}}
Targets = WeightTargetsMirandaOutput(seq_input_file, predicted_targets, scores)
print('computed weighted targets for all genes', len(Targets))

# get CNV gene status
CNV_status = sort_genes_CNV_status(CNV_file)
print('recorded CNV gene status', len(CNV_status))

# add CNV status
for gene in Targets:
    Targets[gene].append(CNV_status[gene])
    assert len(Targets[gene]) == 4, 'gene in Targets does not have all required values'
print('added gene CNV status to each gene')

# create a list of list with weighted targets for CNV and non-CNV genes
AllData = [[], []]
for gene in Targets:
    if Targets[gene][-1] == 'CNV':
        # add weighted number of targets to cnv list
        AllData[0].append(Targets[gene][2])
    elif Targets[gene][-1] == 'not_CNV':
        # add weighted number of targets to non-cnv list
        AllData[1].append(Targets[gene][2])
print('generated lists of target sites for CNV and non-CNV genes')


# perform stattistical tests between CNV and non-CNV genes
Pval = stats.ranksums(AllData[0], AllData[1])[1]
print('cnv: N = {0}, mean = {1}, non-cnv: N = {2}, mean = {3}, P = {4}'.format(len(AllData[0]),
      np.mean(AllData[0]), len(AllData[1]), np.mean(AllData[1]), Pval))    
print('compared CNV and non-CNV genes')

# create figure
fig = plt.figure(1, figsize = (4, 2.5))
# create subplot in figure
# add a plot to figure (N row, N column, plot N)
ax = fig.add_subplot(1, 1, 1)
# create a list of positions for the box plot    
BoxPositions = [0, 0.4]
# use a boxplot
bp = ax.boxplot(AllData, showmeans = True, showfliers = False, widths = 0.3,
                positions = BoxPositions, patch_artist = True) 

i = 0
# change box, whisker color to black
for box in bp['boxes']:
    # change line color
    box.set(color = 'black')
    if i % 2 == 0:
        # CNV data
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
   
# set font for all text in figure
FigFont = {'fontname':'Arial'}   
    
# write label for y and x axis
ax.set_ylabel('Weighted number of miRNA sites per nucleotide', color = 'black',  size = 8, ha = 'center', **FigFont)
ax.set_xlabel('miRNA expression level', color = 'black',  size = 8, ha = 'center', **FigFont)

# create list of labels and tick positions for the X axis
xtickpos = [0.2, 1.1]
CNVLabels = ['CNV', 'non-CNV']
# write label for x axis
plt.xticks(xtickpos, CNVLabels, ha = 'center', fontsize = 8, **FigFont)

## add a range for the Y and X axis
plt.ylim([0, 0.14])
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
   
# annotate Graph with significance level
if Pval >= 0.05:
    Significance = ''
elif Pval < 0.05 and Pval >= 0.01:
    Significance = '*'
elif Pval < 0.01 and Pval >= 0.001:
    Significance = '**'
elif Pval < 0.001:
    Significance = '***'
    
# annotate figure with significance levels
Ypos = 0.13
Xpos = 0.2
ax.text(Xpos, Ypos, Significance, horizontalalignment = 'center', verticalalignment = 'center', color = 'black', size = 8)

## add legend relative to ax1 using ax1 coordinates
#C = mpatches.Patch(facecolor = '#a6bddb', edgecolor = 'black', linewidth = 1, label= 'CNV')
#N = mpatches.Patch(facecolor = '#99d8c9', edgecolor = 'black', linewidth = 1, label= 'non-CNV')
#ax.legend(handles = [C, N], loc = (0.2, 1), fontsize = 8, frameon = False, ncol = 2)

## build outputfile with arguments
#outputfile = 'truc_' + domain + '_' + chromos + '_' + cnv_length
#print(outputfile)

# save figure
fig.savefig('truc.pdf', bbox_inches = 'tight')

