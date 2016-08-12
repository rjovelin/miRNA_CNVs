# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 12:11:59 2016

@author: RJovelin
"""


# use this script to compare target sites between human CNV and non-CNv genes
# for mirnas in different expression groups

# usage PlotTargetSitesExpressionGroup.py [options]
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


# create a dictionary with {species: {mirna accession : mature names pairs}}
AccessionNames = MatchmiRNAAccessionNumbers('name', miRBaseFile = 'miRNA.dat')
# check that mirna names are for human mirna accession numbers
for accession in AccessionNames['Homo_sapiens']:
    for name in AccessionNames['Homo_sapiens'][accession]:
        assert name.startswith('hsa'), 'mirna name should start with hsa'
print('matches accessions with names')

# create a dictionary {mirna accession: expression level}
miRNAEpression = miRBAsemiRNAExpression('mirna_read_count.txt')
print('obtained mirna expression')

# sort the mirna accession numbers according to the mirna expression 
AccessionExpressionGroups = SortmiRNAQuartileExpression('Homo_sapiens', miRNAExpression, AccessionNames)    
print('sorted accessions according to expression')

missing = []
# Sort mirna names in expression groups
LowExp, ModerateExp, MediumExp, HighExp = [], [], [], []
for accession in AccessionNames['Homo_sapiens']:
    if accession in AccessionExpressionGroups[0]:
        # add names to low expression group
        for name in AccessionNames['Homo_sapiens'][accession]:
            LowExp.append(name)
    elif accession in AccessionExpressionGroups[1]:
        # add names to moderate expression group
        for name in AccessionNames['Homo_sapiens'][accession]:
            ModerateExp.append(name)
    elif accession in AccessionExpressionGroups[2]:
        # add names to medium expression group
        for name in AccessionNames['Homo_sapiens'][accession]:
            MediumExp.append(name)
    elif accession in AccessionExpressionGroups[3]:
        # add names to high expression group
        for name in AccessionNames['Homo_sapiens'][accession]:
            HighExp.append(name)
    else:
        missing.append(name)
assert len(missing) == 0, 'some mirnas are missing'


# get the seq input file
seq_input_file = 'H_sapiens_' + domain + '_' + chromos + '_targetscan.txt'
# get the predicted targets output file
predicted_targets = 'H_sapiens_' + domain + '_' + chromos + '_predicted_sites_' + predictors[i] + '.txt'
# use DGV 2015 release 
CNV_file = 'H_sapiens_GRCh37_2015_CNV_all_length_valid_chromos.txt'


# record the number of miranda target sites for each gene for each expression group
#  {gene: [N_targets, Sequence_length, N_targets_normalized, CNV_status}}
TargetsLowExp = SelectmiRNAsMirandaOutput(seq_input_file, predicted_targets, LowExp, 'all')
TargetsModerateExp = SelectmiRNAsMirandaOutput(seq_input_file, predicted_targets, ModerateExp, 'all')
TargetsMediumExp = SelectmiRNAsMirandaOutput(seq_input_file, predicted_targets, MediumExp, 'all')
TargetsHighExp = SelectmiRNAsMirandaOutput(seq_input_file, predicted_targets, HighExp, 'all')
print('recorded targets for each mirnas in each expression group')

# get CNV gene status
CNV_status = sort_genes_CNV_status(CNV_file)
print('recorded CNV gene status')


# add CNV status
for gene in TargetsLowExp:
    TargetsLowExp.append(CNV_status[gene])
    assert len(TargetsLowExp[gene]) == 4, 'gene in TargetsLowExp does not have all required values'
for gene in TargetsModerateExp:
    TargetsModerateExp[gene].append(CNV_status[gene])
    assert len(TargetsModerateExp[gene]) == 4, 'gene in TargetsModerateExp does not have all required values'
for gene in TargetsMediumExp:
    TargetsMediumExp[gene].append(CNV_status[gene])
    assert len(TargetsMediumExp[gene]) == 4, 'gene in TargetsMediumExp does not have all required values'
for gene in TargetsHighExp:
    TargetsHighExp[gene].append(CNV_status[gene])
    assert len(TargetsHighExp[gene]) == 4, 'gene in TargetsHighExp does not have all required values'
print('added gene CNV status to each gene')

# create a dict with expression group as key and a list with the number of targets
# for CNv and non-CNV genes {expression: [[CNV], [non-CNV]]}
TargetsData = []
# initialize dict
TargetsData['low'], TargetsData['moderate'], TargetsData['medium'], TargetsData['high'] = [[], []], [[], []], [[], []], [[], []]
for gene in TargetsLowExp:
    if TargetsLowExp[gene][-1] == 'CNV':
        # add number of normalized targets to cnv list 
        TargetsData['low'][0].append(TargetsLowExp[gene][2])
    elif TargetsLowExp[gene][-1] == 'not_CNV':
        # add number of normalized targets to non-cnv list
        TargetsData['low'][1].append(TargetsLowExp[gene][2])
for gene in TargetsModerateExp:
    if TargetsModerateExp[gene][-1] == 'CNV':
        # add number of normalized targets to cnv list
        TargetsData['moderate'][0].append(TargetsModerateExp[gene][2])
    elif TargetsModerateExp[gene][-1] == 'not_CNV':
        # add number of normalized targets to non-cnv list
        TargetsData['moderate'][1].append(TargetsModerateExp[gene][2])
for gene in TargetsMediumExp:
    if TargetsMediumExp[gene][-1] == 'CNV':
        # add number of normalized targets to cnv list
        TargetsData['medium'][0].append(TargetsMediumExp[gene][2])
    elif TargetsMediumExp[gene][-1] == 'not_CNV':
        # add numner of normnalized targets to non-cnv list
        TargetsData['medium'][1].append(TargetsMediumExp[gene][2])
for gene in TargetsHighExp:
    if TargetsHighExp[gene][-1] == 'CNV':
        # add number of normalized targets to cnv list
        TargetsData['high'][0].append(TargetsHighExp[gene][2])
    elif TargetsHighExp[gene][-1] == 'not_CNV':
        # add number of normalized targets to non-cnv list
        TargetsData['high'][1].append(TargetsHighExp[gene][2])
print('generated lists of target sites for CNV and non-CNV genes')


# perform stattistical tests between CNV and non-CNV genes
# create dicts to store results {expression group: P-value normalized targets}
CompTargets = {}
for group in TargetsData:
    Pval = stats.ranksums(TargetsData[group][0], TargetsData[group][1])[1]
    CompTargets[group] = Pval
print('compared CNV and non-CNV genes')
# print P-values
for group in CompTargets:
    print(group, CompTargets[group])


# make a list with all data
AllData = []
Groups = ['low', 'moderate', 'medium', 'high']
# loop over expression groups,  
for group in Groups:
    # add list of targets for cnv genes
    AllData.append(TargetsData[group][0])
    # add list of targets for non-cnv genes
    AllData.append(TargetsData[group][1])
print('data consolidated in array')


# create figure
fig = plt.figure(1, figsize = (4, 2.5))


# CreateAx(Columns, Rows, Position, Data, figure, Title, YLabel, YMax, SpeciesNames, XScale):

# create subplot in figure
# add a plot to figure (N row, N column, plot N)
ax = figure.add_subplot(1, 1, 1)
# create a list of positions for the box plot    
BoxPositions = [0, 0.4, 0.9, 1.3, 1.8, 2.2, 2.7, 3.1]
# use a boxplot
bp = ax.boxplot(AllData, showmeans = True, showfliers = False, widths = 0.3,
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
   

## write title   
#ax.set_title(Title, size = 8)
    
# set font for all text in figure
FigFont = {'fontname':'Arial'}   
    
# write label for y and x axis
ax.set_ylabel('Normalized number of miRNA\nsites per gene', color = 'black',  size = 8, ha = 'center', **FigFont)
ax.set_xlabel('miRNA expression level', color = 'black',  size = 8, ha = 'center', **FigFont)

# create list of labels and tick positions for the X axis
xtickpos = [0.2, 1.1, 2, 2.9]
ExGroups = ['Low', 'Moderate', 'Medium', 'High']
# write label for x axis
plt.xticks(xtickpos, ExGroups, ha = 'center', fontsize = 8, **FigFont)

## add a range for the Y and X axis
#plt.ylim([0, YMax])
plt.xlim([-0.25, 3.15])
    
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
   
      
## annotate Graph with significance level
#PvalTargetScanAbsolute, PvalMirandaAbsolute, PvalTargetScanNormalized, PvalMirandaNormalized = [], [], [], []
#for species in species_names:
#    # get the significance level for target sites
#    if CompTargetscan[species][0] >= 0.05:
#        PvalTargetScanAbsolute.append('')
#    elif CompTargetscan[species][0] < 0.05 and CompTargetscan[species][0] >= 0.01:
#        PvalTargetScanAbsolute.append('*')
#    elif CompTargetscan[species][0] < 0.01 and CompTargetscan[species][0] >= 0.001:
#        PvalTargetScanAbsolute.append('**')
#    elif CompTargetscan[species][0] < 0.001:
#        PvalTargetScanAbsolute.append('***')
#    # get the significance level of normalized target sites
#    if CompTargetscan[species][1] >= 0.05:
#        PvalTargetScanNormalized.append('')
#    elif CompTargetscan[species][1] < 0.05 and CompTargetscan[species][1] >= 0.01:
#        PvalTargetScanNormalized.append('*')
#    elif CompTargetscan[species][1] < 0.01 and CompTargetscan[species][1] >= 0.001:
#        PvalTargetScanNormalized.append('**')
#    elif CompTargetscan[species][1] < 0.001:
#        PvalTargetScanNormalized.append('***')
#
#
## create list of Y and X positions to annotate figure with significance level
#if domain == '3UTR':
#    # make a list of Y positions
#    YposTargetscanAbsolute = [1390, 420, 400, 990, 390, 580]
#    YposMirandaAbsolute = [820, 270, 210, 590, 220, 320]
#    YposTargetscanNormalized = [0.42, 0.11, 0.16, 0.33, 0.15, 0.21]
#    YposMirandaNormalized = [0.34, 0.10, 0.12, 0.23, 0.13, 0.16]
#    Xpos = [0.2, 1.1, 2, 2.9, 3.8, 4.7]
#
#
#for i in range(len(PvalTargetScanAbsolute)):
#    ax1.text(Xpos[i], YposTargetscanAbsolute[i], PvalTargetScanAbsolute[i], horizontalalignment = 'center',
#             verticalalignment = 'center', color = 'black', size = 8)
#    ax2.text(Xpos[i], YposTargetscanNormalized[i], PvalTargetScanNormalized[i], horizontalalignment = 'center',
#             verticalalignment = 'center', color = 'black', size = 8)         
#for i in range(len(PvalMirandaAbsolute)):
#    ax3.text(Xpos[i], YposMirandaAbsolute[i], PvalMirandaAbsolute[i], horizontalalignment = 'center',
#             verticalalignment = 'center', color = 'black', size = 8)
#    ax4.text(Xpos[i], YposMirandaNormalized[i], PvalMirandaNormalized[i], horizontalalignment = 'center',
#             verticalalignment = 'center', color = 'black', size = 8)
#
## add legend relative to ax1 using ax1 coordinates
#C = mpatches.Patch(facecolor = '#a6cee3', edgecolor = 'black', linewidth = 1, label= 'CNV')
#N = mpatches.Patch(facecolor = '#b2df8a', edgecolor = 'black', linewidth = 1, label= 'non-CNV')
#ax1.legend(handles = [C, N], loc = (0.8, 1.2), fontsize = 8, frameon = False, ncol = 2)




#
## build outputfile with arguments
#outputfile = 'truc_' + domain + '_' + chromos + '_' + cnv_length
#print(outputfile)

# save figure
fig.savefig('truc.pdf', bbox_inches = 'tight')

