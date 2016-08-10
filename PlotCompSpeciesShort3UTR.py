# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 11:59:52 2016

@author: RJovelin
"""


# use this script to compare the number of miRNA targets per nucleotide
# between human and each other vertebrate species for genes with short 3'UTRs
# plot a box plot for target sites in CDS and 5'UTR for CNV genes and for non-CNV genes

# usage PlotCompSpeciesShort3UTR.py

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
# import custom modules
from CNV_miRNAs import *

# keep genes on assembled nuclear chromosomes
chromos = 'valid_chromos'
keep_valid_chromos = True
# consider all CNVs
cnv_length = 'CNV_all_length'
CNV_size = 'all'
# use only mirnas that are shared between species pairs
conservation = 'shared_miRs'
# consider short 3'UTR is L < 7bp
L = 7
# alternatives
# chromos = 'all_chromos'
# cnv_length = 'CNV_greater_1Kb'
# CNV_size = 'long'
# conservation = all_miRs


# get the number of target sites for CNV genes and for non-CNV genes in human and each other vertebrates

# use the 2015 release of the DGV
human_CNV_file = 'H_sapiens_GRCh37_2015_CNV_all_length_valid_chromos.txt'

# make a dict with species name : species code pairs
species_codes = {'H_sapiens': 'Hsa', 'P_troglodytes': 'Ptr', 'M_mulatta': 'Mml',
                 'M_musculus': 'Mmu', 'B_taurus': 'Bta', 'G_gallus':'Gga'}
# make a list of species names
SpeciesNames = ['H_sapiens', 'P_troglodytes', 'M_mulatta', 'M_musculus', 'B_taurus', 'G_gallus']

# create a dict for human-species pairwise comparison of target sites for CNV genes or non-CNV genes
# and for sites in CDS or 5'UTR of genes with short 3'UTRs
# {species_name : {domain: {'CNV_status': [[sites_Hsa], [sites_sp2]]}}}
HsaSpeciesTargets = {}


# loop over domains
for domain in ['5UTR', 'CDS']:
    # loop over non-human species
    for i in range(1, len(SpeciesNames)):
        # compare human to each other species
        print(SpeciesNames[0], SpeciesNames[i])
        # get genome files
        genome_sp1 = SpeciesNames[0] + '_genome.txt'
        genome_sp2 = SpeciesNames[i] + '_genome.txt'
        # get chromos files
        valid_chromos_sp1 = SpeciesNames[0] + '_valid_chromos.txt'
        valid_chromos_sp2 = SpeciesNames[i] + '_valid_chromos.txt'
        # get GFF annotation files
        GFF_sp1 = SpeciesNames[0] + '.gff3'
        GFF_sp2 = SpeciesNames[i] + '.gff3'
    
        # get targetscan sequence input file
        seq_input_sp1 = SpeciesNames[0] + '_' + domain + '_' + chromos + '_targetscan.txt'
        seq_input_sp2 = SpeciesNames[i] + '_' + domain + '_' + chromos + '_targetscan.txt'
           
        # get predictor output
        predicted_targets_sp1 = SpeciesNames[0] + '_' + domain + '_' + chromos + '_predicted_sites_targetscan.txt'
        predicted_targets_sp2 = SpeciesNames[i] + '_' + domain + '_' + chromos + '_predicted_sites_targetscan.txt'
                    
        # get UTR files
        UTR_sp1 = SpeciesNames[0] + '_3UTR_length_' + chromos + '.txt'      
        UTR_sp2 = SpeciesNames[i] + '_3UTR_length_' + chromos + '.txt'
             
        # get CNV files
        CNV_file_sp1 = human_CNV_file
        CNV_file_sp2 = SpeciesNames[i] + '_' + cnv_length + '_' + chromos + '.txt'
                
        # get mature mirna files
        mature_sp1 = SpeciesNames[0] + '_mature.txt'
        mature_sp2 = SpeciesNames[i] + '_mature.txt'        
            
        # get a set of shared seeds between species pair
        shared_seeds = find_conserved_mirna_families(mature_sp1, mature_sp2)
        
        # get CNV gene status
        CNV_status1 = sort_genes_CNV_status(CNV_file_sp1)
        CNV_status2 = sort_genes_CNV_status(CNV_file_sp2)
           
        # parse predictor outputfile {gene: [N_sites, seq_length, N_sites/seq_length]}
        # consider only shared mirna families between species pairs
        # filter targets for shared mirnas
        targets_sp1 = parse_targetscan_output(seq_input_sp1, predicted_targets_sp1, 'conserved', shared_seeds, mature_sp1)
        targets_sp2 = parse_targetscan_output(seq_input_sp2, predicted_targets_sp2, 'conserved', shared_seeds, mature_sp2)
            
        # sort genes by UTR length
        UTR_length_sp1 = sort_genes_3UTR_length(UTR_sp1, L)
        UTR_length_sp2 = sort_genes_3UTR_length(UTR_sp2, L)
            
        # create list of normalized targets for CNV and non-CNV genes
        # record only genes with short 3' UTR and sort genes that are CNV or not_CNV    
        Sp1TargetsCNV, Sp2TargetsCNV, Sp1TargetsNonCNV, Sp2TargetsNonCNV =  [], [], [], []
        # record only genes with short 3' UTR and that are CNV    
        Sp1TargetsCNV = [targets_sp1[gene][2] for gene in targets_sp1 if UTR_length_sp1[gene] == 'short' and CNV_status1[gene] == 'CNV']    
        Sp2TargetsCNV = [targets_sp2[gene][2] for gene in targets_sp2 if UTR_length_sp2[gene] == 'short' and CNV_status2[gene] == 'CNV']   
        # record only genes with short 3' UTR and that not_CNV    
        Sp1TargetsNonCNV = [targets_sp1[gene][2] for gene in targets_sp1 if UTR_length_sp1[gene] == 'long' and CNV_status1[gene] == 'not_CNV']
        Sp2TargetsNonCNV = [targets_sp2[gene][2] for gene in targets_sp2 if UTR_length_sp2[gene] == 'long' and CNV_status2[gene] == 'not_CNV']    
        print('got mirna targets for short CNV and non-CNV genes')
        print('{0} and {1} CNV genes for {2} and {3}'.format(len(Sp1TargetsCNV), len(Sp2TargetsCNV), SpeciesNames[0], SpeciesNames[1]))
    
        # populate dict
        HsaSpeciesTargets[SpeciesNames[i]] = {}
        HsaSpeciesTargets[SpeciesNames[i]][domain] = {}
        # copy lists to avoid modifying list values in dict
        HsaSpeciesTargets[SpeciesNames[i]][domain]['CNV'] = [Sp1TargetsCNV[:], Sp2TargetsCNV[:]]
        HsaSpeciesTargets[SpeciesNames[i]][domain]['not_CNV'] = [Sp1TargetsNonCNV[:], Sp2TargetsNonCNV[:]]
    


# plot the number of target sites of mirna families shared between human and
# with each other species for CNV and for non-CNV genes
 
# create paralell lists for target sites in CNV and non-CNV genes located in 5'UTR or CDS 
CNVGenesCDS, CNVGenes5UTR, NonCNVGenesCDS, NonCNVGenes5UTR = [], [], [], [] 
# loop over species names
for species in SpeciesNames[1:]:
    # loop over domain
    for domain in HsaSpeciesTargets[species]:
        # loop over CNV status
        for status in HsaSpeciesTargets[species][domain]:
            # check domain and CNV status
            if domain == 'CDS' and status == 'CNV':
                for i in range(2):
                    CNVGenesCDS.append(HsaSpeciesTargets[species][domain][status][i])
            elif domain == 'CDS' and status =='not_CNV':
                for i in range(2):
                    NonCNVGenesCDS.append(HsaSpeciesTargets[species][domain][status][i])
            elif domain == '5UTR' and status == 'CNV':
                for i in range(2):
                    CNVGenes5UTR.append(HsaSpeciesTargets[species][domain][status][i])
            elif domain == '5UTR' and status == 'not_CNV':
                for i in range(2):
                    NonCNVGenes5UTR.append(HsaSpeciesTargets[species][domain][status][i])
print('generated lists of target sites for CNV and non-CNV genes')


## perform stattistical tests between human and other species for CNV and non-CNV genes
## create dicts to store results {species: [P-value absolute targets, P-value normalized targets]}
#CompTargetsCDS, CompTargets5UTR = {}, {}
#for species in HsaSpeciesTargets:
#    for domain in HsaSpeciesTargets[species]:
#        
#
#
#
#
#
#
# SpeciesDataTargetscanAbsolute:
#    Pabs = stats.ranksums(SpeciesDataTargetscanAbsolute[species][0], SpeciesDataTargetscanAbsolute[species][1])[1]
#    Pnorm = stats.ranksums(SpeciesDataTargetscanNormalized[species][0], SpeciesDataTargetscanNormalized[species][1])[1]
#    CompTargetscan[species] = [Pabs, Pnorm]
#for species in SpeciesDataMirandaAbsolute:
#    Pabs = stats.ranksums(SpeciesDataMirandaAbsolute[species][0], SpeciesDataMirandaAbsolute[species][1])[1]    
#    Pnorm = stats.ranksums(SpeciesDataMirandaNormalized[species][0], SpeciesDataMirandaNormalized[species][1])[1]    
#    CompMiranda[species] = [Pabs, Pnorm]
#print('compared CNV and non-CNV genes')
#
## print P-values
#for species in CompTargetscan:
#    print('targetscan', species, CompTargetscan[species])
#for species in CompMiranda:
#    print('miranda', species, CompMiranda[species])



# create figure
fig = plt.figure(1, figsize = (8, 5))

# create list of labels and tick positions for the X axis
xtickpos = [0.2, 1.1, 2, 2.9, 3.8]
Names = [species_codes[i] for i in SpeciesNames]
print(Names)

# create a function to format the subplots
def CreateAx(Columns, Rows, Position, Data, figure, Title, SpeciesNames, XScale):
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
    BoxPositions = [0, 0.4, 0.9, 1.3, 1.8, 2.2, 2.7, 3.1, 3.6, 4]
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
    plt.ylim([0, 0.16])
    
    plt.xlim([-0.25, 4.25])
    
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



# create a lit of species abbrviations
labelnames = [species_codes[i] for i in SpeciesNames[1:]]

# plot box plots
ax1 = CreateAx(2, 2, 1, CNVGenes5UTR, fig, '5\'UTRs of CNV genes', labelnames, xtickpos)
ax2 = CreateAx(2, 2, 2, NonCNVGenes5UTR, fig, '5\'UTRs of non-CNV genes', labelnames, xtickpos)
ax3 = CreateAx(2, 2, 3, CNVGenesCDS, fig, 'CDS of CNV genes', labelnames, xtickpos)
ax4 = CreateAx(2, 2, 4, NonCNVGenesCDS, fig, 'CDS of non-CNV genes', labelnames, xtickpos)

# add subplot labels
ax1.text(-1.3, 0.18, 'A', horizontalalignment = 'center', verticalalignment = 'center', color = 'black', size = 10)
ax1.text(5.5, 0.18, 'B', horizontalalignment = 'center', verticalalignment = 'center', color = 'black', size = 10)   
ax3.text(-1.3, 0.18, 'C', horizontalalignment = 'center', verticalalignment = 'center', color = 'black', size = 10)
ax3.text(5.5, 0.18, 'D', horizontalalignment = 'center', verticalalignment = 'center', color = 'black', size = 10) 


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
#for species in species_names:
#    # get the significance level for target sites
#    if CompMiranda[species][0] >= 0.05:
#        PvalMirandaAbsolute.append('')
#    elif CompMiranda[species][0] < 0.05 and CompMiranda[species][0] >= 0.01:
#        PvalMirandaAbsolute.append('*')
#    elif CompMiranda[species][0] < 0.01 and CompMiranda[species][0] >= 0.001:
#        PvalMirandaAbsolute.append('**')
#    elif CompMiranda[species][0] < 0.001:
#        PvalMirandaAbsolute.append('***')
#    # get the significance level of normalized target sites
#    if CompMiranda[species][1] >= 0.05:
#        PvalMirandaNormalized.append('')
#    elif CompMiranda[species][1] < 0.05 and CompMiranda[species][1] >= 0.01:
#        PvalMirandaNormalized.append('*')
#    elif CompMiranda[species][1] < 0.01 and CompMiranda[species][1] >= 0.001:
#        PvalMirandaNormalized.append('**')
#    elif CompMiranda[species][1] < 0.001:
#        PvalMirandaNormalized.append('***')
#
## create list of Y and X positions to annotate figure with significance level
#if domain == '3UTR':
#    # make a list of Y positions
#    YposTargetscanAbsolute = [1390, 420, 400, 990, 390, 580]
#    YposMirandaAbsolute = [820, 270, 210, 590, 220, 320]
#    YposTargetscanNormalized = [0.42, 0.11, 0.16, 0.33, 0.15, 0.21]
#    YposMirandaNormalized = [0.34, 0.10, 0.12, 0.23, 0.13, 0.16]
#    Xpos = [0.2, 1.1, 2, 2.9, 3.8, 4.7]
#elif domain == 'CDS':
#    # make a list of Y positions
#    YposTargetscanAbsolute = [1400, 420, 400, 1000, 390, 590]
#    YposMirandaAbsolute = [820, 280, 210, 600, 220, 320]
#    YposTargetscanNormalized = [0.42, 0.11, 0.16, 0.33, 0.16, 0.21]
#    YposMirandaNormalized = [0.35, 0.10, 0.12, 0.24, 0.13, 0.16]
#    Xpos = [0.2, 1.1, 2, 2.9, 3.8, 4.7]
#elif domain == '5UTR':
#    # make a list of Y positions
#    YposTargetscanAbsolute = [4000, 1490, 2500, 2020, 1020, 610]
#    YposMirandaAbsolute = [2520, 900, 1500, 1500, 300, 400]
#    YposTargetscanNormalized = [0.40, 0.12, 0.16, 0.33, 0.15, 0.23]
#    YposMirandaNormalized = [0.33, 0.09, 0.12, 0.23, 0.12, 0.17]
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


# make sure subplots do not overlap
plt.tight_layout()


## build outputfile with arguments
#outputfile = 'truc_' + domain + '_' + chromos + '_' + cnv_length
#print(outputfile)

# save figure
fig.savefig('truc.pdf', bbox_inches = 'tight')
