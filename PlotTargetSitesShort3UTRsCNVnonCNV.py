# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 16:12:06 2016

@author: RJovelin
"""


# use this script to compare the number target sites per nucleotide between CNV and non-CNV genes
# that have short 3'UTRs in each species 

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


# usage python3 PlotTargetSitesShort3UTRsCNVnonCNV.py options
# [5UTR/CDS]: use target sites in 5'UTR or CDS of short 3'UTR genes
# [7/15]: 3'UTR length, genes with 3'UTR length < 7bp or < 15bp are considered short 3'UTR genes

# get the region to consider to predict target sites [5UTr or CDS]
domain = sys.argv[1]
print(domain)
# get the minimum 3'UTR length
L = int(sys.argv[2])
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

# create a function to format the subplots
def CreateAx(Columns, Rows, Position, Data, figure, Title, YMax, SpeciesNames, XScale):
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
    # create a list of positions for the box plot    
    BoxPositions = [0, 0.4, 0.9, 1.3, 1.8, 2.2, 2.7, 3.1, 3.6, 4, 4.5, 4.9]
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
    
    
    # write title   
    ax.set_title(Title, size = 8)
    
    # set font for all text in figure
    FigFont = {'fontname':'Arial'}   
    
    # write label for y axis
    ax.set_ylabel('Normalized number of miRNA\nsites per gene', color = 'black',  size = 8, ha = 'center', **FigFont)

    # write label for x axis
    plt.xticks(XScale, SpeciesNames, ha = 'center', fontsize = 8, **FigFont)

    # add a range for the Y axis
    plt.ylim([0, YMax])
    
    plt.xlim([-0.25, 5.15])
    
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


# plot data for targetscan
ax1 = CreateAx(2, 1, 1, AllDataTargetscan, fig, 'TargetScan', 0.45, Names, xtickpos)
ax2 = CreateAx(2, 1, 2, AllDataMiranda, fig, 'miRanda', 0.45, Names, xtickpos)

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


# create list of Y and X positions to annotate figure with significance level
if domain == 'CDS':
    # make a list of Y positions
    YposTargetscan = [0.41, 0.11, 0.16, 0.28, 0.14, 0.21]
    YposMiranda = [0.32, 0.08, 0.11, 0.19, 0.11, 0.16]
    Xpos = [0.2, 1.1, 2, 2.9, 3.8, 4.7]
elif domain == '5UTR':
    # make a list of Y positions
    YposTargetscan = [0.41, 0.12, 0.17, 0.30, 0.16, 0.24]
    YposMiranda = [0.32, 0.09, 0.12, 0.21, 0.12, 0.17]
    Xpos = [0.2, 1.1, 2, 2.9, 3.8, 4.7]

for i in range(len(PvalTargetScan)):
    ax1.text(Xpos[i], YposTargetscan[i], PvalTargetScan[i], horizontalalignment = 'center',
             verticalalignment = 'center', color = 'black', size = 8)
for i in range(len(PvalMiranda)):
    ax2.text(Xpos[i], YposMiranda[i], PvalMiranda[i], horizontalalignment = 'center',
             verticalalignment = 'center', color = 'black', size = 8)

# make sure subplots do not overlap
plt.tight_layout()

# get outputfile
if L == 15:
    outputfile = 'PlotShort3UTR_15bp_' + domain + '_' + chromos + '_' + cnv_length + '_NormalizedSites' 
elif L == 7:
    outputfile = 'PlotShort3UTR_7bp_' + domain + '_' + chromos + '_' + cnv_length + '_NormalizedSites' 
print(outputfile)

# save figure
fig.savefig(outputfile + '.eps', bbox_inches = 'tight')
       
