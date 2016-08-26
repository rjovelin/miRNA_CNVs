# -*- coding: utf-8 -*-
"""
Created on Fri Aug 19 12:11:53 2016

@author: RJovelin
"""

# use this script to compare target sites between CNV and non-CNv genes
# for mirnas in different expression groups in each species

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


# make a list of species names
SpeciesNames = ['H_sapiens', 'M_mulatta', 'M_musculus', 'B_taurus', 'G_gallus']
# create a dict with full species names
Genus = {'H_sapiens': 'Homo_sapiens', 'M_mulatta': 'Macaca_mulatta', 'M_musculus': 'Mus_musculus',
         'B_taurus': 'Bos_taurus', 'G_gallus':'Gallus_gallus'}
# make a dictionary of species names : species code
species_codes = {'H_sapiens': 'Hsa', 'M_mulatta': 'Mml', 'M_musculus': 'Mmu',
                 'B_taurus': 'Bta', 'G_gallus':'Gga'}

# create a dictionary with {species: {mirna accession : mature names pairs}}
AccessionNames = MatchmiRNAAccessionNumbers('name', miRBaseFile = 'miRNA.dat')
# check that mirna names: mirna accession numbers pairs correspond to the correct species
for species in SpeciesNames:
    print(Genus[species])
    for accession in AccessionNames[Genus[species]]:
        for name in AccessionNames[Genus[species]][accession]:
            if Genus[species] == 'Homo_sapiens':
                assert name.startswith('hsa'), 'mirna name should start with hsa'
            elif Genus[species] == 'Macaca_mulatta':
                assert name.startswith('mml'), 'mirna name should start with mml'
            elif Genus[species] == 'Mus_musculus':
                assert name.startswith('mmu'), 'mirna name should start with mmu'
            elif Genus[species] == 'Gallus_gallus':
                assert name.startswith('gga'), 'mirna name should start with gga'
            elif Genus[species] == 'Bos_taurus':
                assert name.startswith('bta'), 'mirna name should start with bta'
print('matches accessions with names')

# create a dictionary {mirna accession: expression level}
miRNAExpression = miRBAsemiRNAExpression('mirna_read_count.txt')
print('obtained mirna expression')

# note that chimp has no expressionm, so chimp is removed from analyses
ChimpExpression = [miRNAExpression[mirna] for mirna in miRNAExpression if mirna in AccessionNames['Pan_troglodytes']]
print('no expression recorded in miRBase for Chimp: ', len(ChimpExpression))


# create a dict with expression group as key and a list with the number of targets
# for CNv and non-CNV genes for all species {species: {expression: [[CNV], [non-CNV]]}}
TargetsData = {}
# loop over species, sort mirna to expression group and record the number of targets for each groups and CNV status
for species in SpeciesNames:
    # get the genus_species
    Species = Genus[species]
    print(Species)
    # sort the mirna accession numbers according to the mirna expression 
    AccessionExpressionGroups = SortmiRNAQuartileExpression(Species, miRNAExpression, AccessionNames)  
    print('sorted accessions according to expression')
    # Sort mirna names in expression groups
    missing = set()
    LowExp, ModerateExp, MediumExp, HighExp = [], [], [], []
    for accession in AccessionNames[Species]:
        if accession in AccessionExpressionGroups[0]:
            # add names to low expression group
            for name in AccessionNames[Species][accession]:
                LowExp.append(name)
        elif accession in AccessionExpressionGroups[1]:
            # add names to moderate expression group
            for name in AccessionNames[Species][accession]:
                ModerateExp.append(name)
        elif accession in AccessionExpressionGroups[2]:
            # add names to medium expression group
            for name in AccessionNames[Species][accession]:
                MediumExp.append(name)
        elif accession in AccessionExpressionGroups[3]:
            # add names to high expression group
            for name in AccessionNames[Species][accession]:
                HighExp.append(name)
        else:
            missing.add(accession)
    print('{0} miRNAs without expression in {1}'.format(len(missing), Species))

    # get the seq input file
    seq_input_file = species + '_' + domain + '_' + chromos + '_targetscan.txt'
    # get the predicted targets output file
    predicted_targets = species + '_' + domain + '_' + chromos + '_predicted_sites_miranda.txt'
    # get the CNV file, use DGV 2015 release for human
    if species == 'H_sapiens':
        CNV_file = 'H_sapiens_GRCh37_2015_CNV_all_length_valid_chromos.txt'
    else:
        CNV_file = species + '_' + cnv_length + '_' + chromos + '.txt'         
        
    # record the number of miranda target sites for each gene for each expression group
    #  {gene: [N_targets, Sequence_length, N_targets_normalized, CNV_status}}
    TargetsLowExp = SelectmiRNAsMirandaOutput(seq_input_file, predicted_targets, LowExp, 'all')
    print('recorded targets for each mirnas in low expression group')
    TargetsModerateExp = SelectmiRNAsMirandaOutput(seq_input_file, predicted_targets, ModerateExp, 'all')
    print('recorded targets for each mirnas in moderate expression group')
    TargetsMediumExp = SelectmiRNAsMirandaOutput(seq_input_file, predicted_targets, MediumExp, 'all')
    print('recorded targets for each mirnas in medium expression group')
    TargetsHighExp = SelectmiRNAsMirandaOutput(seq_input_file, predicted_targets, HighExp, 'all')
    print('recorded targets for each mirnas in high expression group')

    # get CNV gene status
    CNV_status = sort_genes_CNV_status(CNV_file)
    print('recorded CNV gene status')

    # add CNV status
    for gene in TargetsLowExp:
        TargetsLowExp[gene].append(CNV_status[gene])
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

    # initialize dict
    TargetsData[species] = {}
    TargetsData[species]['low'], TargetsData[species]['moderate'], TargetsData[species]['medium'], TargetsData[species]['high'] = [[], []], [[], []], [[], []], [[], []]
    for gene in TargetsLowExp:
        if TargetsLowExp[gene][-1] == 'CNV':
            # add number of normalized targets to cnv list 
            TargetsData[species]['low'][0].append(TargetsLowExp[gene][2])
        elif TargetsLowExp[gene][-1] == 'not_CNV':
            # add number of normalized targets to non-cnv list
            TargetsData[species]['low'][1].append(TargetsLowExp[gene][2])
    for gene in TargetsModerateExp:
        if TargetsModerateExp[gene][-1] == 'CNV':
            # add number of normalized targets to cnv list
            TargetsData[species]['moderate'][0].append(TargetsModerateExp[gene][2])
        elif TargetsModerateExp[gene][-1] == 'not_CNV':
            # add number of normalized targets to non-cnv list
            TargetsData[species]['moderate'][1].append(TargetsModerateExp[gene][2])
    for gene in TargetsMediumExp:
        if TargetsMediumExp[gene][-1] == 'CNV':
            # add number of normalized targets to cnv list
            TargetsData[species]['medium'][0].append(TargetsMediumExp[gene][2])
        elif TargetsMediumExp[gene][-1] == 'not_CNV':
            # add numner of normnalized targets to non-cnv list
            TargetsData[species]['medium'][1].append(TargetsMediumExp[gene][2])
    for gene in TargetsHighExp:
        if TargetsHighExp[gene][-1] == 'CNV':
            # add number of normalized targets to cnv list
            TargetsData[species]['high'][0].append(TargetsHighExp[gene][2])
        elif TargetsHighExp[gene][-1] == 'not_CNV':
            # add number of normalized targets to non-cnv list
            TargetsData[species]['high'][1].append(TargetsHighExp[gene][2])
    print('sorted targets for CNV and non-CNV genes in each expression group for {0}'.format(species))


# perform stattistical tests between CNV and non-CNV genes
# create dicts to store results {species: {expression group: P-value normalized targets}}
CompTargets = {}
for species in TargetsData:
    CompTargets[species] = {}
    for group in TargetsData[species]:
        Pval = stats.ranksums(TargetsData[species][group][0], TargetsData[species][group][1])[1]
        CompTargets[species][group] = Pval    
print('compared CNV and non-CNV genes')


# make a list of expression group 
Groups = ['low', 'moderate', 'medium', 'high']

# {species: [[targets CNV low], [targets non-CNV high], [ targets CNV moderate], [targets non-CNV moderate],
#            [targets CNV medium], [targets non-CNV medium], [targets CNV high], [targets non-CNV high]}
AllData = {}
for species in TargetsData:
    AllData[species] = []
    for group in Groups:
        AllData[species].append(TargetsData[species][group][0])
        AllData[species].append(TargetsData[species][group][1])
print('data consolidated in array for each species')



# create list of labels and tick positions for the X axis
#xtickpos = [0.2, 1.1, 2, 2.9, 3.8, 4.7]

# create a function to format the subplots
def CreateAx(Columns, Rows, Position, Data, figure, Title, YMax, YAxisLine):
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
    BoxPositions = [0, 0.4, 0.9, 1.3, 1.8, 2.2, 2.7, 3.1]

    # use a boxplot
    bp = ax.boxplot(Data, showmeans = True, showfliers = False, widths = 0.3,
                    positions = BoxPositions, patch_artist = True) 

    # color CNV and non-CNV boxes differently
    CNVColor = ['#a6bddb','#74a9cf','#2b8cbe','#045a8d']
    NonCNVColor = ['#99d8c9','#66c2a4','#2ca25f','#006d2c']
    i, j = 0, 0    
    # change box, whisker color to black
    for box in bp['boxes']:
        # change line color
        box.set(color = 'black')
        if i % 2 == 0:
            # CNV data
            box.set(facecolor = CNVColor[j])
        else:
            box.set(facecolor = NonCNVColor[j])
        i += 1
        j = int(i / 2)
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
    # write label for x and y axis
    ax.set_ylabel('Normalized number of miRNA\nsites per gene', color = 'black',  size = 8, ha = 'center', **FigFont)
    if YAxisLine == True:
        ax.set_xlabel('miRNA expression level', color = 'black',  size = 8, ha = 'center', **FigFont)
    
    # write label for x axis
    plt.xticks([0.2, 1.1, 2, 2.9], list(map(lambda x: x.capitalize(), Groups)), ha = 'center', fontsize = 8, **FigFont)
    # add a range for the Y axis
    plt.ylim([0, YMax])
    plt.xlim([-0.25, 3.35])

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


# create figure
fig = plt.figure(1, figsize = (3.5, 7))

# plot data, note that chimp has no expression data 
ax1 = CreateAx(1, 5, 1, AllData[SpeciesNames[0]], fig, Genus[SpeciesNames[0]].replace('_', ' '), 0.5, False)
ax2 = CreateAx(1, 5, 2, AllData[SpeciesNames[1]], fig, Genus[SpeciesNames[1]].replace('_', ' '), 0.5, False)
ax3 = CreateAx(1, 5, 3, AllData[SpeciesNames[2]], fig, Genus[SpeciesNames[2]].replace('_', ' '), 0.5, False)
ax4 = CreateAx(1, 5, 4, AllData[SpeciesNames[3]], fig, Genus[SpeciesNames[3]].replace('_', ' '), 0.5, False)
ax5 = CreateAx(1, 5, 5, AllData[SpeciesNames[4]], fig, Genus[SpeciesNames[4]].replace('_', ' '), 0.5, True)

# annotate Graph with significance level
Pvalues = {}
for species in SpeciesNames:
    Pvalues[species] = []
    for group in Groups:
        # get the significance level for target sites
        if CompTargets[species][group] >= 0.05:
            Pvalues[species].append('')
        elif CompTargets[species][group] < 0.05 and CompTargets[species][group] >= 0.01:
            Pvalues[species].append('*')
        elif CompTargets[species][group] < 0.01 and CompTargets[species][group] >= 0.001:
            Pvalues[species].append('**')
        elif CompTargets[species][group] < 0.001:
            Pvalues[species].append('***')

# create list of Y and X positions to annotate figure with significance level
if domain == '3UTR':
    # make a list of Y positions
    Ypos = [0.13, 0.10, 0.09, 0.09]
    Xpos = [0.2, 1.1, 2, 2.9]
elif domain == '5UTR':
    # make a list of Y positions
    Ypos = [0.13, 0.10, 0.09, 0.09]
    Xpos = [0.2, 1.1, 2, 2.9]
elif domain == 'CDS':
    # make a list of Y positions
    Ypos = [0.13, 0.10, 0.09, 0.09]
    Xpos = [0.2, 1.1, 2, 2.9]


## annotate figure with significance levels
#for i in range(len(Pvalues)):
#    ax.text(Xpos[i], Ypos[i], Pvalues[i], horizontalalignment = 'center', verticalalignment = 'center', color = 'black', size = 8)


# add legend relative to ax1 using ax1 coordinates
C = mpatches.Patch(facecolor = '#a6bddb', edgecolor = 'black', linewidth = 1, label= 'CNV')
N = mpatches.Patch(facecolor = '#99d8c9', edgecolor = 'black', linewidth = 1, label= 'non-CNV')
ax.legend(handles = [C, N], loc = (0.2, 1.2), fontsize = 8, frameon = False, ncol = 2)






## build outputfile with arguments
#outputfile = 'PlotTargetsmiRNAExpression_' + domain + '_' + chromos + '_' + cnv_length
#print(outputfile)

# save figure
fig.savefig('truc.pdf', bbox_inches = 'tight')

