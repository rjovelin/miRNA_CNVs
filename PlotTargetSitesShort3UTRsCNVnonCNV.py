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
        UTR_length = {}
        infile = open(UTR_file)
        infile.readline()
        for line in infile:
            line = line.rstrip()
            if line != '':
                line = line.split('\t')
                # get gene name, and 3' UTR length
                gene, L3UTR = line[0], int(line[2])
                # populate dict
                if L3UTR < L:
                    UTR_length[gene] = 'short'
                elif L3UTR > L:
                    UTR_length[gene] = 'long'
        infile.close()                    
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
    for gene in targetscan:
        if targetscan[species][gene][-1] == 'CNV':
            SpeciesDataTargetscan[species][0].append(targetscan[species][gene][2])
        elif targetscan[species][gene][-1] == 'not_CNV':
            SpeciesDataTargetscan[species][1].append(targetscan[species][gene][2])
for species in miranda:
    # initialize list values
    SpeciesDataMiranda[species] = [[], []]
    # populate inner lists with number of mirna target sites per nucleotide
    for gene in miranda:
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
fig = plt.figure(1, figsize = (5, 2.5))

    
# create a function to format the subplots
def CreateAx(Columns, Rows, Position, Xscale, Data, figure, Title, XLabel = False, YLabel = False):
    # create subplot in figure
    # add a plot to figure (N row, N column, plot N)
    ax = figure.add_subplot(Rows, Columns, Position)
    # create a list of positions for the box plot    
    BoxPositions = [0, 0.6, 0.8, 1.4, 1.6, 2.2, 2.4, 3, 3.2, 3.8, 4, 4.6]
    # use a boxplot
    bp = ax.boxplot(Data, showmeans = True, showfliers = False, widths = 0.5,
                    positions = BoxPositions, patch_artist = True) 
    
    
          
    
    
    
    
    
    
    
    # color CNV boxes in grey
    i = 0    
    # change box, whisker color to black
    for box in bp['boxes']:
        # change line color
        box.set(color = 'black')
        if i % 2 == 0:
            # CNV data, color box in grey
            box.set(facecolor = "grey")
        else:
            box.set(facecolor = 'white')
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
        mean.set(marker = 'o', markeredgecolor = 'black', markerfacecolor = 'black', markersize = 4)
    
    
    
    
    
    ax.set_title(Title, size = 10)
    
    
    
    
    
    # write label for y axis
    ytext = ax.set_ylabel('Normalized number of miRNA\nsites per gene', color = 'grey',  size = 10, ha = 'center', fontname = 'Arial', family = 'sans-serif')










# set tick label
names = []
for species in species_names:
    names.append(species_codes[species] + '_CNV')
    names.append(species_codes[species] + '_non_CNV')

# add labels to x-ticks, rotate and align right, set size to 14
ax.set_xticklabels(names, rotation = 30, ha = 'right', size = 10, fontname = 'Arial', family = 'sans-serif')


# add a light grey horizontal grid to the plot, semi-transparent, 
ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)  

# hide these grids behind plot objects
ax.set_axisbelow(True)



    
    

# do not show lines around figure  
ax.spines["top"].set_visible(False)    
ax.spines["bottom"].set_visible(False)    
ax.spines["right"].set_visible(False)    
ax.spines["left"].set_visible(False)  

# create a list with range of x-axis values
xvals = [i + 0.5 for i in range(len(names) + 1)]
# Set a buffer around the edge of the x-axis
plt.xlim([min(xvals)- 0.5, max(xvals)+ 0.5])

# do not show ticks
plt.tick_params(
    axis='y',       # changes apply to the x-axis and y-axis (other option : x, y)
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    right = 'off',
    left = 'off',          
    labelbottom='off', # labels along the bottom edge are off 
    colors = 'grey',
    labelsize = 10)
      

# do not show ticks
plt.tick_params(
    axis='x',       # changes apply to the x-axis and y-axis (other option : x, y)
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    right = 'off',
    left = 'off',          
    labelbottom='on', # labels along the bottom edge are on 
    colors = 'black',
    labelsize = 10)

# Set the tick labels font name
for label in ax.get_yticklabels():
    label.set_fontname('Arial')

# add title
if predictor == 'targetscan':
    plt.title('TargetScan', size = 10, fontname = 'Arial')  
elif predictor == 'miranda':
    plt.title('miRanda', size = 10, fontname = 'Arial')


# get outputfile
outputfile = 'Fig_short3UTR_CNVvsNonCNV_' + domain + '_' + chromos + '_' + cnv_length + '_' + predictor + '_normalized_sites' 
print(outputfile)

# save figure
fig.savefig(outputfile + '.eps', bbox_inches = 'tight')
       
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    if YLabel == True:
        # set y axis label
        ax.set_ylabel('Protein pairs', size = 7, ha = 'center', fontname = 'Helvetica', family = 'sans-serif')
    if XLabel == True:
        # set x axis label
        ax.set_xlabel('Protein distance', size = 7, ha = 'center', fontname = 'Helvetica', family = 'sans-serif')

    # do not show lines around figure, keep bottow line  
    ax.spines["top"].set_visible(False)    
    ax.spines["right"].set_visible(False)    
    ax.spines["left"].set_visible(False)      
    
    if BottomLine == True:
        ax.spines["bottom"].set_visible(True)
    else:
        ax.spines["bottom"].set_visible(False)
    
        # add a light grey horizontal grid to the plot, semi-transparent, 
    ax.yaxis.grid(True, linestyle=':', which='major', color='lightgrey', alpha=0.5, linewidth = 0.5)  
    # hide these grids behind plot objects
    ax.set_axisbelow(True)    
        
    if BottomLine == True:
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
            labelsize = 7,
            direction = 'out') # ticks are outside the frame when bottom = 'on'  
    else:
        # do not show ticks
        plt.tick_params(
            axis='both',       # changes apply to the x-axis and y-axis (other option : x, y)
            which='both',      # both major and minor ticks are affected
            bottom='off',      # ticks along the bottom edge are off
            top='off',         # ticks along the top edge are off
            right = 'off',
            left = 'off',          
            labelbottom='off', # labels along the bottom edge are on
            colors = 'black',
            labelsize = 7,
            direction = 'out') # ticks are outside the frame when bottom = 'on'  

    if XLabel == True:
        # determine tick position on x axis
        xpos =  [i /10  for i in range(11)]
        Dist = ['0', '', '', '', '', '0.5', '', '', '', '', '1']
        # set up tick positions and labels
        plt.xticks(xpos, Dist, rotation = 0, fontsize = 7, ha = 'center', fontname = 'Helvetica')

    if YLabel == True:
        plt.yticks(fontsize = 8)
  
    return ax

















# set width for all subplots    
width = 0.1

ax1 = CreateAx(4, 5, 1, [i / 10 for i in range(10)], Distances[families[0]], fig, families[0], width, YLabel = True)
ax2 = CreateAx(4, 5, 2, [i / 10 for i in range(10)], Distances[families[1]], fig, families[1], width)
ax3 = CreateAx(4, 5, 3, [i / 10 for i in range(10)], Distances[families[2]], fig, families[2], width)
ax4 = CreateAx(4, 5, 4, [i / 10 for i in range(10)], Distances[families[3]], fig, families[3], width)
ax5 = CreateAx(4, 5, 5, [i / 10 for i in range(10)], Distances[families[4]], fig, families[4], width, YLabel = True)
ax6 = CreateAx(4, 5, 6, [i / 10 for i in range(10)], Distances[families[5]], fig, families[5], width)
ax7 = CreateAx(4, 5, 7, [i / 10 for i in range(10)], Distances[families[6]], fig, families[6], width)
ax8 = CreateAx(4, 5, 8, [i / 10 for i in range(10)], Distances[families[7]], fig, families[7], width)
ax9 = CreateAx(4, 5, 9, [i / 10 for i in range(10)], Distances[families[8]], fig, families[8], width, YLabel = True)
ax10 = CreateAx(4, 5, 10, [i / 10 for i in range(10)], Distances[families[9]], fig, families[9], width)
ax11 = CreateAx(4, 5, 11, [i / 10 for i in range(10)], Distances[families[10]], fig, families[10], width)
ax12 = CreateAx(4, 5, 12, [i / 10 for i in range(10)], Distances[families[11]], fig, families[11], width)
ax13 = CreateAx(4, 5, 13, [i / 10 for i in range(10)], Distances[families[12]], fig, families[12], width, YLabel = True)
ax14 = CreateAx(4, 5, 14, [i / 10 for i in range(10)], Distances[families[13]], fig, families[13], width)
ax15 = CreateAx(4, 5, 15, [i / 10 for i in range(10)], Distances[families[14]], fig, families[14], width)
ax16 = CreateAx(4, 5, 16, [i / 10 for i in range(10)], Distances[families[15]], fig, families[15], width)
ax17 = CreateAx(4, 5, 17, [i / 10 for i in range(10)], Distances[families[16]], fig, families[16], width, XLabel = True, YLabel = True, BottomLine = True)
ax18 = CreateAx(4, 5, 18, [i / 10 for i in range(10)], Distances[families[17]], fig, families[17], width, XLabel = True, BottomLine = True)
ax19 = CreateAx(4, 5, 19, [i / 10 for i in range(10)], Distances[families[18]], fig, families[18], width, XLabel = True, BottomLine = True)

# make sure subplots do not overlap
plt.tight_layout()

# save figure
fig.savefig('ChemoProtDivergenceFamilies.pdf', bbox_inches = 'tight')


#################





# write label for y axis
ytext = ax.set_ylabel('Normalized number of miRNA\nsites per gene', color = 'grey',  size = 10, ha = 'center', fontname = 'Arial', family = 'sans-serif')

# set tick label
names = []
for species in species_names:
    names.append(species_codes[species] + '_CNV')
    names.append(species_codes[species] + '_non_CNV')

# add labels to x-ticks, rotate and align right, set size to 14
ax.set_xticklabels(names, rotation = 30, ha = 'right', size = 10, fontname = 'Arial', family = 'sans-serif')


# add a light grey horizontal grid to the plot, semi-transparent, 
ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)  

# hide these grids behind plot objects
ax.set_axisbelow(True)

# remove top axes and right axes ticks
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()


# use a boxplot
bp = ax.boxplot(all_data, showmeans = True, showfliers = False, widths = 0.7, labels = names, patch_artist = True) 
    
# color CNV boxes in grey
i = 0    
# change box, whisker color to black
for box in bp['boxes']:
    # change line color
    box.set(color = 'black')
    if i % 2 == 0:
        # CNV data, color box in grey
        box.set(facecolor = "grey")
    else:
        box.set(facecolor = 'white')
    i += 1
        
# change whisker color ro black
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
    mean.set(marker = 'o', markeredgecolor = 'black', markerfacecolor = 'black', markersize = 4)
    

# do not show lines around figure  
ax.spines["top"].set_visible(False)    
ax.spines["bottom"].set_visible(False)    
ax.spines["right"].set_visible(False)    
ax.spines["left"].set_visible(False)  

# create a list with range of x-axis values
xvals = [i + 0.5 for i in range(len(names) + 1)]
# Set a buffer around the edge of the x-axis
plt.xlim([min(xvals)- 0.5, max(xvals)+ 0.5])

# do not show ticks
plt.tick_params(
    axis='y',       # changes apply to the x-axis and y-axis (other option : x, y)
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    right = 'off',
    left = 'off',          
    labelbottom='off', # labels along the bottom edge are off 
    colors = 'grey',
    labelsize = 10)
      

# do not show ticks
plt.tick_params(
    axis='x',       # changes apply to the x-axis and y-axis (other option : x, y)
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    right = 'off',
    left = 'off',          
    labelbottom='on', # labels along the bottom edge are on 
    colors = 'black',
    labelsize = 10)

# Set the tick labels font name
for label in ax.get_yticklabels():
    label.set_fontname('Arial')

# add title
if predictor == 'targetscan':
    plt.title('TargetScan', size = 10, fontname = 'Arial')  
elif predictor == 'miranda':
    plt.title('miRanda', size = 10, fontname = 'Arial')


# get outputfile
outputfile = 'Fig_short3UTR_CNVvsNonCNV_' + domain + '_' + chromos + '_' + cnv_length + '_' + predictor + '_normalized_sites' 
print(outputfile)

# save figure
fig.savefig(outputfile + '.eps', bbox_inches = 'tight')
    
    
    
    
    
    
   