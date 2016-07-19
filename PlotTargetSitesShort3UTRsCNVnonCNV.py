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

#######################

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
        

################# continue here        
        
        
        
        
        
        
        
        
        
        
        
    
        # get CNV gene status
        CNV_status = sort_genes_CNV_status(CNV_file)
        print('CNV status', len(CNV_status))
    
        # write number of target sites, sequence length, normalized targets and CNV status to file
        newfile = open(outputfile, 'w')
        # write header
        newfile.write('\t'.join(['Gene', 'N_targets', 'Sequence_length', 'N_targets_normalized', 'CNV_status']) + '\n')
    
        # loop over genes in targets
        for gene in targets:
            # record only genes with short 3'UTR
            if gene in UTR_length and UTR_length[gene] == 'short':
                # write number of sites, sequence length and normalized number of targets to file
                newfile.write('\t'.join([gene, str(targets[gene][0]), str(targets[gene][1]), str(targets[gene][2])]) + '\t')
                # write CNV status
                newfile.write(CNV_status[gene] + '\n')
        # close file after writing
        newfile.close()
        print('done writing targets to file')



















##########################


# create a dict of {species name : [[CNV], [non-CNV]]}
species_data = {}

# loop over species names
# parse the summary target file, extract number of normalized sites and CNV status
# populate dict, use dict to generate figure        
for species in species_names:
    print(species)
    
    # get the summary target file
    if species == 'H_sapiens':
        # use DGV 2015 release 
        summary_file = species + '_short3UTRgenes_' + domain + '_summary_' + predictor + '_' + chromos + '_' + cnv_length + '_GRCh37_2015.txt'
    else:
        summary_file = species + '_short3UTRgenes_' + domain + '_summary_' + predictor + '_' + chromos + '_' + cnv_length + '.txt'
    
    # parse the summary table into a list of list of values
    # [CNV_targets, nonCNV_targets, CNV_seq, nonCNV_seq, CNV_normalized, nonCNV_normalized]    
    parsed_data = parse_summary_table_targets(summary_file)
    
    # populate dicts with normalized target sites
    species_data[species] = [parsed_data[4], parsed_data[5]]
    print(species, len(species_data[species][0]), len(species_data[species][1]))
        
print('nb species', len(species_data))    
    
# make a list of data
all_data = []
# loop over species name, append list of CNV values and list of non-CNV values
for species in species_names:
    all_data.append(species_data[species][0])
    all_data.append(species_data[species][1])     
print('data consolidated in array')
    


# create figure
fig = plt.figure(1, figsize = (4.3,2.56))

# add a plot to figure (1 row, 1 column, 1 plot)
ax = fig.add_subplot(1, 1, 1)    

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
    
    
    
    
    
    
   