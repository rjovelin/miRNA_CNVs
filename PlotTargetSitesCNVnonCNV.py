# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 15:16:57 2016

@author: RJovelin
"""


# use this script to compare target sites and normalized target sites between CNV and non-CNv genes in each species

# usage PlotTargetSitesCNVnonCNV.py [options]
# [3UTR/5UTR/CDS]: gene domain to consider
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



# get the region to consider to predict target sites [3UTR or 5UTr or CDS]
domain = sys.argv[1]
print(domain)
# check the format of the figure (eps by default or pdf if specified)
if len(sys.argv) == 2:
    extension = '.eps'
elif len(sys.argv) == 3:
    assert sys.argv[2] == 'pdf'
    extension = '.pdf'

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
            
        # initialize inner dict
        if i == 0:
            targetscan[species] = {}
            for gene in targets:
                targetscan[species][gene] = list(targets[gene])
        elif i == 1:
            miranda[species] = {}
            for gene in targets:
                miranda[species][gene] = list(targets[gene])
        
# check that the same number of species are recorded for miranda and targetscan
assert len(targetscan) == len(miranda), 'different number of species depending on predictor'
# check that list values have the correct number of items
for species in targetscan:
    for gene in targetscan[species]:
        assert len(targetscan[species][gene]) == 4, 'gene in targetscan does not have all required values'
    for gene  in miranda[species]:
        assert len(miranda[species][gene]) == 4, 'gene in miranda does not have all required values'


# create a dicts with species as keys and list of targets for CNV and non-CNV genes,
# with absolute number of targets and with normalized number of targets {species name : [[CNV], [non-CNV]]}
SpeciesDataTargetscanAbsolute, SpeciesDataMirandaAbsolute  = {}, {}
SpeciesDataTargetscanNormalized, SpeciesDataMirandaNormalized = {}, {}

# loop over species names, get the number of sites for CNV and non-CNV genes
for species in targetscan:
    # initialise list value
    SpeciesDataTargetscanAbsolute[species] = [[], []]
    SpeciesDataTargetscanNormalized[species] = [[], []]
    # populate inner lists with number of miRNA target sites per nucleotide
    for gene in targetscan[species]:
        if targetscan[species][gene][-1] == 'CNV':
            SpeciesDataTargetscanAbsolute[species][0].append(targetscan[species][gene][0])
            SpeciesDataTargetscanNormalized[species][0].append(targetscan[species][gene][2])
        elif targetscan[species][gene][-1] == 'not_CNV':
            SpeciesDataTargetscanAbsolute[species][1].append(targetscan[species][gene][0])
            SpeciesDataTargetscanNormalized[species][1].append(targetscan[species][gene][2])
for species in miranda:
    # initialize list values
    SpeciesDataMirandaAbsolute[species] = [[], []]
    SpeciesDataMirandaNormalized[species] = [[], []]
    # populate inner lists with number of mirna target sites per nucleotide
    for gene in miranda[species]:
        if miranda[species][gene][-1] == 'CNV':
            SpeciesDataMirandaAbsolute[species][0].append(miranda[species][gene][0])
            SpeciesDataMirandaNormalized[species][0].append(miranda[species][gene][2])
        elif miranda[species][gene][-1] == 'not_CNV':
            SpeciesDataMirandaAbsolute[species][1].append(miranda[species][gene][0])
            SpeciesDataMirandaNormalized[species][1].append(miranda[species][gene][2])
print('generated lists of target sites for CNV and non-CNV genes')


# perform stattistical tests between CNV and non-CNV genes
# create dicts to store results {species: [P-value absolute targets, P-value normalized targets]}
CompTargetscan, CompMiranda = {}, {}
for species in SpeciesDataTargetscanAbsolute:
    Pabs = stats.ranksums(SpeciesDataTargetscanAbsolute[species][0], SpeciesDataTargetscanAbsolute[species][1])[1]
    Pnorm = stats.ranksums(SpeciesDataTargetscanNormalized[species][0], SpeciesDataTargetscanNormalized[species][1])[1]
    CompTargetscan[species] = [Pabs, Pnorm]
for species in SpeciesDataMirandaAbsolute:
    Pabs = stats.ranksums(SpeciesDataMirandaAbsolute[species][0], SpeciesDataMirandaAbsolute[species][1])[1]    
    Pnorm = stats.ranksums(SpeciesDataMirandaNormalized[species][0], SpeciesDataMirandaNormalized[species][1])[1]    
    CompMiranda[species] = [Pabs, Pnorm]
print('compared CNV and non-CNV genes')

# print P-values
for species in CompTargetscan:
    print('targetscan', species, CompTargetscan[species])
for species in CompMiranda:
    print('miranda', species, CompMiranda[species])


# make a list of data for each predictor
AllDataTargetscanAbsolute, AllDataMirandaAbsolute, AllDataTargetscanNormalized, AllDataMirandaNormalized = [], [], [], []

# make a list of species names to loop from
species_names = ['H_sapiens', 'P_troglodytes', 'M_mulatta', 'M_musculus', 'B_taurus', 'G_gallus']
# loop over species in species names list and populate data lists, keeping the same order for targetscan and miranda
for species in species_names:
    # append list of target sites for CNV genes
    AllDataTargetscanAbsolute.append(SpeciesDataTargetscanAbsolute[species][0])
    AllDataTargetscanNormalized.append(SpeciesDataTargetscanNormalized[species][0])
    AllDataMirandaAbsolute.append(SpeciesDataMirandaAbsolute[species][0])
    AllDataMirandaNormalized.append(SpeciesDataMirandaNormalized[species][0])    
    # append list of target sites for non-CNV genes
    AllDataTargetscanAbsolute.append(SpeciesDataTargetscanAbsolute[species][1])
    AllDataTargetscanNormalized.append(SpeciesDataTargetscanNormalized[species][1])
    AllDataMirandaAbsolute.append(SpeciesDataMirandaAbsolute[species][1])
    AllDataMirandaNormalized.append(SpeciesDataMirandaNormalized[species][1])
print('data consolidated in array')


# create figure
fig = plt.figure(1, figsize = (8, 5))

# create list of labels and tick positions for the X axis
#xtickpos = [0.35, 1.25, 2.15, 3.05, 3.95, 4.85]
xtickpos = [0.2, 1.1, 2, 2.9, 3.8, 4.7]
Names = [species_codes[i] for i in species_names]
print(Names)

# create a function to format the subplots
def CreateAx(Columns, Rows, Position, Data, figure, Title, YLabel, YMax, SpeciesNames, XScale):
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
    ax.set_ylabel(YLabel, color = 'black',  size = 8, ha = 'center', **FigFont)

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
if domain == '5UTR':
    ax1 = CreateAx(2, 2, 1, AllDataTargetscanAbsolute, fig, 'TargetScan', 'Number of miRNA sites per gene', 4000, Names, xtickpos)
    ax2 = CreateAx(2, 2, 2, AllDataTargetscanNormalized, fig, 'TargetScan', 'Normalized number of miRNA\nsites per gene', 0.45, Names, xtickpos)
    ax3 = CreateAx(2, 2, 3, AllDataMirandaAbsolute, fig, 'miRanda', 'Number of miRNA sites per gene', 2600, Names, xtickpos)
    ax4 = CreateAx(2, 2, 4, AllDataMirandaNormalized, fig, 'miRanda', 'Normalized number of miRNA\nsites per gene', 0.40,  Names, xtickpos)
else:
    ax1 = CreateAx(2, 2, 1, AllDataTargetscanAbsolute, fig, 'TargetScan', 'Number of miRNA sites per gene', 1400, Names, xtickpos)
    ax2 = CreateAx(2, 2, 2, AllDataTargetscanNormalized, fig, 'TargetScan', 'Normalized number of miRNA\nsites per gene', 0.45, Names, xtickpos)
    ax3 = CreateAx(2, 2, 3, AllDataMirandaAbsolute, fig, 'miRanda', 'Number of miRNA sites per gene', 800, Names, xtickpos)
    ax4 = CreateAx(2, 2, 4, AllDataMirandaNormalized, fig, 'miRanda', 'Normalized number of miRNA\nsites per gene', 0.45,  Names, xtickpos)


# add subplot label
if domain == '5UTR':
    ax1.text(-1.3, 1550, 'A', horizontalalignment = 'center', verticalalignment = 'center', color = 'black', size = 10)
    ax1.text(5.5, 1550, 'B', horizontalalignment = 'center', verticalalignment = 'center', color = 'black', size = 10)   
    ax3.text(-1.3, 2700, 'C', horizontalalignment = 'center', verticalalignment = 'center', color = 'black', size = 10)
    ax3.text(5.5, 2700, 'D', horizontalalignment = 'center', verticalalignment = 'center', color = 'black', size = 10) 
else:
    ax1.text(-1.3, 1550, 'A', horizontalalignment = 'center', verticalalignment = 'center', color = 'black', size = 10)
    ax1.text(5.5, 1550, 'B', horizontalalignment = 'center', verticalalignment = 'center', color = 'black', size = 10)   
    ax3.text(-1.3, 880, 'C', horizontalalignment = 'center', verticalalignment = 'center', color = 'black', size = 10)
    ax3.text(5.5, 880, 'D', horizontalalignment = 'center', verticalalignment = 'center', color = 'black', size = 10) 


# annotate Graph with significance level
PvalTargetScanAbsolute, PvalMirandaAbsolute, PvalTargetScanNormalized, PvalMirandaNormalized = [], [], [], []
for species in species_names:
    # get the significance level for target sites
    if CompTargetscan[species][0] >= 0.05:
        PvalTargetScanAbsolute.append('')
    elif CompTargetscan[species][0] < 0.05 and CompTargetscan[species][0] >= 0.01:
        PvalTargetScanAbsolute.append('*')
    elif CompTargetscan[species][0] < 0.01 and CompTargetscan[species][0] >= 0.001:
        PvalTargetScanAbsolute.append('**')
    elif CompTargetscan[species][0] < 0.001:
        PvalTargetScanAbsolute.append('***')
    # get the significance level of normalized target sites
    if CompTargetscan[species][1] >= 0.05:
        PvalTargetScanNormalized.append('')
    elif CompTargetscan[species][1] < 0.05 and CompTargetscan[species][1] >= 0.01:
        PvalTargetScanNormalized.append('*')
    elif CompTargetscan[species][1] < 0.01 and CompTargetscan[species][1] >= 0.001:
        PvalTargetScanNormalized.append('**')
    elif CompTargetscan[species][1] < 0.001:
        PvalTargetScanNormalized.append('***')

for species in species_names:
    # get the significance level for target sites
    if CompMiranda[species][0] >= 0.05:
        PvalMirandaAbsolute.append('')
    elif CompMiranda[species][0] < 0.05 and CompMiranda[species][0] >= 0.01:
        PvalMirandaAbsolute.append('*')
    elif CompMiranda[species][0] < 0.01 and CompMiranda[species][0] >= 0.001:
        PvalMirandaAbsolute.append('**')
    elif CompMiranda[species][0] < 0.001:
        PvalMirandaAbsolute.append('***')
    # get the significance level of normalized target sites
    if CompMiranda[species][1] >= 0.05:
        PvalMirandaNormalized.append('')
    elif CompMiranda[species][1] < 0.05 and CompMiranda[species][1] >= 0.01:
        PvalMirandaNormalized.append('*')
    elif CompMiranda[species][1] < 0.01 and CompMiranda[species][1] >= 0.001:
        PvalMirandaNormalized.append('**')
    elif CompMiranda[species][1] < 0.001:
        PvalMirandaNormalized.append('***')

# create list of Y and X positions to annotate figure with significance level
if domain == '3UTR':
    # make a list of Y positions
    YposTargetscanAbsolute = [1390, 420, 400, 990, 390, 580]
    YposMirandaAbsolute = [820, 270, 210, 590, 220, 320]
    YposTargetscanNormalized = [0.42, 0.11, 0.16, 0.33, 0.15, 0.21]
    YposMirandaNormalized = [0.34, 0.10, 0.12, 0.23, 0.13, 0.16]
    Xpos = [0.2, 1.1, 2, 2.9, 3.8, 4.7]
elif domain == 'CDS':
    # make a list of Y positions
    YposTargetscanAbsolute = [1400, 420, 400, 1000, 390, 590]
    YposMirandaAbsolute = [820, 280, 210, 600, 220, 320]
    YposTargetscanNormalized = [0.42, 0.11, 0.16, 0.33, 0.16, 0.21]
    YposMirandaNormalized = [0.35, 0.10, 0.12, 0.24, 0.13, 0.16]
    Xpos = [0.2, 1.1, 2, 2.9, 3.8, 4.7]
elif domain == '5UTR':
    # make a list of Y positions
    YposTargetscanAbsolute = [4000, 1490, 2500, 2020, 1020, 610]
    YposMirandaAbsolute = [2520, 900, 1500, 1500, 300, 400]
    YposTargetscanNormalized = [0.40, 0.12, 0.16, 0.33, 0.15, 0.23]
    YposMirandaNormalized = [0.33, 0.09, 0.12, 0.23, 0.12, 0.17]
    Xpos = [0.2, 1.1, 2, 2.9, 3.8, 4.7]


for i in range(len(PvalTargetScanAbsolute)):
    ax1.text(Xpos[i], YposTargetscanAbsolute[i], PvalTargetScanAbsolute[i], horizontalalignment = 'center',
             verticalalignment = 'center', color = 'black', size = 8)
    ax2.text(Xpos[i], YposTargetscanNormalized[i], PvalTargetScanNormalized[i], horizontalalignment = 'center',
             verticalalignment = 'center', color = 'black', size = 8)         
for i in range(len(PvalMirandaAbsolute)):
    ax3.text(Xpos[i], YposMirandaAbsolute[i], PvalMirandaAbsolute[i], horizontalalignment = 'center',
             verticalalignment = 'center', color = 'black', size = 8)
    ax4.text(Xpos[i], YposMirandaNormalized[i], PvalMirandaNormalized[i], horizontalalignment = 'center',
             verticalalignment = 'center', color = 'black', size = 8)

# add legend relative to ax1 using ax1 coordinates
C = mpatches.Patch(facecolor = '#a6cee3', edgecolor = 'black', linewidth = 1, label= 'CNV')
N = mpatches.Patch(facecolor = '#b2df8a', edgecolor = 'black', linewidth = 1, label= 'non-CNV')
ax1.legend(handles = [C, N], loc = (0.8, 1.2), fontsize = 8, frameon = False, ncol = 2)


# make sure subplots do not overlap
plt.tight_layout()


# build outputfile with arguments
outputfile = 'PlotSitesCNVvsNonCNV_' + domain + '_' + chromos + '_' + cnv_length
print(outputfile)

# save figure
fig.savefig(outputfile + extension, bbox_inches = 'tight')

