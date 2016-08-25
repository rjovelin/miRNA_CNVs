# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 15:55:09 2016

@author: RJovelin
"""

# use this script to compare target sites between CNV and non-CNv genes in each species
# assigning different scores to miRNAs according to their expression level

# usage PlotTargetSitesmiRNAScore.py 

# use Agg backend on server without X server
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rc
rc('mathtext', default='regular')
## activate latex text rendering
#rc('text', usetex=True)
# import modules
import numpy as np
from scipy import stats
import math
import os
import sys
# import custom modules
from CNV_miRNAs import *


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


# make a list of species names
SpeciesNames = ['H_sapiens', 'P_troglodytes', 'M_mulatta', 'M_musculus', 'B_taurus', 'G_gallus']
# create a dict with full species names
Genus = {'H_sapiens': 'Homo_sapiens', 'P_troglodytes': 'Pan_troglodytes', 'M_mulatta': 'Macaca_mulatta',
                 'M_musculus': 'Mus_musculus', 'B_taurus': 'Bos_taurus', 'G_gallus':'Gallus_gallus'}
# make a dictionary of species names : species code
species_codes = {'H_sapiens': 'Hsa', 'P_troglodytes': 'Ptr', 'M_mulatta': 'Mml',
                 'M_musculus': 'Mmu', 'B_taurus': 'Bta', 'G_gallus':'Gga'}
# create a list of domains
regions = ['3UTR', '5UTR', 'CDS']


# record weighted number of sites for each species and domain 
# {species: [[3UTR CNV], [3UTR non-CNV], [5UTR CNV], [5UTR non-CNV], [CDS CNV], [CDS non-CNV]]}
AllData = {}
# create a dict to record the scores of each mature mirna for each species {species: {mature name: score}}
for species in SpeciesNames:
    AllData[species] = [[], [], [], [], [], []]
Scores = {}
# loop over species
for species in SpeciesNames:
    print(species)
    # get the file with mature names and accession numbers
    fastafile = species + '_miRBaseMatureAccession.txt'
    # create a dict mature name : score pairs
    scores = TargetScore(fastafile, 1000, Genus[species], miRBaseFile = 'miRNA.dat', ExpressionFile = 'mirna_read_count.txt')    
    print('match scores to mature miRNAs', species, len(scores))
    # populate score dict
    if len(scores) != 0:
        Scores[species] = dict(scores)
        # get the CNV file, use the DGV 2015 release for human
        if species == 'H_sapiens':
            CNV_file = 'H_sapiens_GRCh37_2015_CNV_all_length_valid_chromos.txt'
        else:
            CNV_file = species + '_' + cnv_length + '_' + chromos + '.txt' 
        # get CNV gene status
        CNV_status = sort_genes_CNV_status(CNV_file)
        print('recorded CNV gene status', len(CNV_status))
        # loop over domain
        for domain in regions:
            # get the seq input file
            seq_input_file = species + '_' + domain + '_' + chromos + '_targetscan.txt'
            # get the predicted targets output file
            predicted_targets = species + '_' + domain + '_' + chromos + '_predicted_sites_miranda.txt'
            # record the number of miranda target sites for each gene weighted by the mirna expression score
            #  {gene: [N_targets, Sequence_length, N_targets_normalized, CNV_status}}
            Targets = WeightTargetsMirandaOutput(seq_input_file, predicted_targets, Scores[species])
            print('computed weighted targets for all genes', len(Targets))
            # add CNV status
            for gene in Targets:
                Targets[gene].append(CNV_status[gene])
                assert len(Targets[gene]) == 4, 'gene in Targets does not have all required values'
            print('added gene CNV status to each gene')
            # add the number of weighted targets to CNV and non-CNV lists
            for gene in Targets:
                if domain == '3UTR' and Targets[gene][-1] == 'CNV':
                    # add weigted number of targets to cnv list for 3'UTR
                    AllData[species][0].append(Targets[gene][2])
                elif domain == '3UTR' and Targets[gene][-1] == 'not_CNV':
                    # add weighted number of targets to non-cnv list for 3'UTR
                    AllData[species][1].append(Targets[gene][2])
                elif domain == '5UTR' and Targets[gene][-1] == 'CNV':
                    # add weighted number of targets to cnv for 5'UTR 
                    AllData[species][2].append(Targets[gene][2])
                elif domain == '5UTR' and Targets[gene][-1] == 'not_CNV':
                    # add weighted number of targets to non-cnv for 5'UTR
                    AllData[species][3].append(Targets[gene][2])
                elif domain == 'CDS' and Targets[gene][-1] == 'CNV':
                    # add weighted number of targets to cnv list for CDS
                    AllData[species][4].append(Targets[gene][2])
                elif domain == 'CDS' and Targets[gene][-1] == 'not_CNV':
                    # add weighted number of targets to non-cnv list for CDS
                    AllData[species][5].append(Targets[gene][2])
            print('generated lists of target sites for {0}'.format(species))
print('generated lists of sites for CNV and non-CNV genes for all species')



# create figure
fig = plt.figure(1, figsize = (4, 7.5))

# create list of labels and tick positions for the X axis
xtickpos = [0.2, 1.1, 2, 2.9, 3.8, 4.7]



# create a function to format the subplots
def CreateAx(Columns, Rows, Position, Data, figure, Title, YLabel, YMax, Domains, XScale, YAxisLine):
    '''
    (int, int, int, dict, figure_object, str, str, int, list, list, bool)
    Take the number of a column, and rows in the figure object and the position of
    the ax in figure, a list of data, a title, a maximum value for the Y axis,
    a list with species names and list of X axis tick positions and return an
    ax instance in the figure
    '''    
    
    
    # create subplot in figure
    # add a plot to figure (N row, N column, plot N)
    ax = figure.add_subplot(Rows, Columns, Position)
    # create a list of positions for the box plot    
    BoxPositions = [0, 0.4, 0.9, 1.3, 1.8, 2.2]
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
    ax.set_title(Title, size = 8, style = 'italic')
    # set font for all text in figure
    FigFont = {'fontname':'Arial'}   
    # write label for y axis
    ax.set_ylabel(YLabel, color = 'black',  size = 8, ha = 'center', **FigFont)
    # write label for x axis
    plt.xticks(XScale, Domains, ha = 'center', fontsize = 8, **FigFont)
    # add a range for the Y axis
    plt.ylim([0, YMax])
    plt.xlim([-0.25, 2.45])
    
    # do not show lines around figure  
    ax.spines["top"].set_visible(False)    
    ax.spines["right"].set_visible(False)    
    ax.spines["left"].set_visible(False)  
    if YAxisLine == False:
        ax.spines["bottom"].set_visible(False)    
    elif YAxisLine == True:
        ax.spines["bottom"].set_visible(True)    
        # offset the spines
        for spine in ax.spines.values():
            spine.set_position(('outward', 5))
    
    # add a light grey horizontal grid to the plot, semi-transparent, 
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5, linewidth = 0.5)  
    # hide these grids behind plot objects
    ax.set_axisbelow(True)


    if YAxisLine == True:
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
    elif YAxisLine == False:
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
            labelsize = 8,
            direction = 'out') # ticks are outside the frame when bottom = 'on'  
    
    # Set the tick labels font name
    for label in ax.get_yticklabels():
        label.set_fontname('Arial')
    # create a margin around the x axis
    plt.margins(0.05)
    
    return ax      

Ylabel = 'Weighted number of\nsites / nt'
# plot data, note that chimp has no expression data 
ax1 = CreateAx(1, 5, 1, AllData[SpeciesNames[0]], fig, Genus[SpeciesNames[0]].replace('_', ' '), Ylabel, 12, regions, xtickpos, False)
ax2 = CreateAx(1, 5, 2, AllData[SpeciesNames[2]], fig, Genus[SpeciesNames[2]].replace('_', ' '), Ylabel, 4, regions, xtickpos, False)
ax3 = CreateAx(1, 5, 3, AllData[SpeciesNames[3]], fig, Genus[SpeciesNames[3]].replace('_', ' '), Ylabel, 8, regions, xtickpos, False)
ax4 = CreateAx(1, 5, 4, AllData[SpeciesNames[4]], fig, Genus[SpeciesNames[4]].replace('_', ' '), Ylabel, 4, regions, xtickpos, False)
ax5 = CreateAx(1, 5, 5, AllData[SpeciesNames[5]], fig, Genus[SpeciesNames[5]].replace('_', ' '), Ylabel, 6, regions, xtickpos, True)

# perform statistical tests between CNV and non-CNV genes
Pval = {}
for species in AllData:
    # initialize list value
    Pval[species] = []
    for i in range(0, 6, 2):
        P = stats.ranksums(AllData[species][i], AllData[species][i+1])[1]    
        Pval[species].append(P)
print('compared CNV and non-CNV genes')

# use stars to show significance levels
Significance = {}
for species in Pval:
    Significance[species] = []
    for i in range(len(Pval[species])):
        if Pval[species][i] >= 0.05:
            Significance[species].append('')
        elif Pval[species][i] < 0.05 and Pval[species][i] >= 0.01:
            Significance[species].append('*')
        elif Pval[species][i] < 0.01 and Pval[species][i] >= 0.001:
            Significance[species].append('**')
        elif Pval[species][i] < 0.001:
            Significance[species].append('***')
print('assessed significance for each comparisons')

# annotate figure with significance level
# create list of positions for significance
Xpos = [0.2, 1.1, 2]
for species in Significance:
    if species == 'H_sapiens':
        # make a list of Y axis position
        Ypos = [11, 11, 11]
        for i in range(len(Ypos)):
            ax1.text(Xpos[i], Ypos[i], Significance[species][i], horizontalalignment = 'center', verticalalignment = 'center', color = 'black', size = 8)
    elif species == 'M_mulatta':
        Ypos = [3.5, 3.5, 3.5]
        for i in range(len(Ypos)):
            ax2.text(Xpos[i], Ypos[i], Significance[species][i], horizontalalignment = 'center', verticalalignment = 'center', color = 'black', size = 8)
    elif species == 'M_musculus':
        Ypos = [7, 7, 7]
        for i in range(len(Ypos)):
            ax3.text(Xpos[i], Ypos[i], Significance[species][i], horizontalalignment = 'center', verticalalignment = 'center', color = 'black', size = 8)
    elif species == 'B_tautus':
        Ypos = [3.7, 3.7, 3.5]
        for i in range(len(Ypos)):
            ax4.text(Xpos[i], Ypos[i], Significance[species][i], horizontalalignment = 'center', verticalalignment = 'center', color = 'black', size = 8)
    elif species == 'G_gallus':
        Ypos = [5, 5, 4]
        for i in range(len(Ypos)):
            ax5.text(Xpos[i], Ypos[i], Significance[species][i], horizontalalignment = 'center', verticalalignment = 'center', color = 'black', size = 8)

## add legend relative to ax1 using ax1 coordinates
#C = mpatches.Patch(facecolor = '#a6bddb', edgecolor = 'black', linewidth = 1, label= 'CNV')
#N = mpatches.Patch(facecolor = '#99d8c9', edgecolor = 'black', linewidth = 1, label= 'non-CNV')
#ax.legend(handles = [C, N], loc = (0.2, 1), fontsize = 8, frameon = False, ncol = 2)

## build outputfile with arguments
#outputfile = 'truc_' + domain + '_' + chromos + '_' + cnv_length
#print(outputfile)


# make sure subplots do not overlap
plt.tight_layout()

# save figure
fig.savefig('truc.pdf', bbox_inches = 'tight')

