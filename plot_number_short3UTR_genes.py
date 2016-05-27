# -*- coding: utf-8 -*-
"""
Created on Sun Nov  8 10:18:44 2015

@author: Richard
"""


# plot the number of genes with short 3'UTR in each species or DVG release as a bar graph


# usage python3 plot_number_short3UTR_genes.py options
# [species/human] number of genes for each species (DVG-2015 for human) or for each DGV release
# [True/False] valid_chromos or all chromos
# [long_CNVs/all_CNVs] CNV > 1 Kb or all CNVs



import os
import sys
import matplotlib.pyplot as plt
from matplotlib import patches as mpatches
import numpy as np

# get the option to consider species or DGV
species = sys.argv[1]
print(species)

# get the option to keep genes on all chromos (False) or only on assembled nuclear chromosomes (True) 
keep_valid_chromos = sys.argv[2]
if keep_valid_chromos == 'True':
    chromos = 'valid_chromos'
elif keep_valid_chromos == 'False':
    chromos = 'all_chromos'
print(keep_valid_chromos, chromos)

# get the type of CNVs to consider from the command [long_CNVs or all_CNVs]
long_CNV = sys.argv[3]
if long_CNV == 'all_CNVs':
    cnv_length = 'CNV_all_length'
elif long_CNV == 'long_CNVs':
    cnv_length = 'CNV_greater_1Kb'
print(long_CNV, cnv_length)

# get the input file with gene counts
if species == 'species':
    counts_file = 'Gene_Counts_Short_3UTR_' + cnv_length + '_' + chromos + '.txt'
elif species == 'human':
    counts_file = 'Gene_Counts_DGV_release_short_3UTR_' + cnv_length + '_' + chromos + '.txt'
print(counts_file)


# open file for reading
infile = open(counts_file, 'r')
# skip 2 first lines
infile.readline()
infile.readline()

# create a dict to store the {species/DVG_release :[total, cnv, non-cnv]} gene counts
gene_counts = {}
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split()
        gene_counts[line[0]] = [int(line[1]), int(line[2]), int(line[3])]
# close file after reading
infile.close()


# make a list of species names or release versions
if species == 'species':
    names = ['H_sapiens', 'P_troglodytes', 'M_mulatta', 'M_musculus', 'B_taurus', 'G_gallus']
    labelnames = ['Hsa', 'Ptr', 'Mml', 'Mmu', 'Bta', 'Gga']
elif species == 'human':
    names = ['GRCh37_2013-05', 'GRCh37_2013-07', 'GRCh37_2014', 'GRCh37_2015']
    labelnames = ['2013a', '2013b', '2014', '2015']
    


# create list of counts for cnv genes parallel to labelnames list
cnv_genes = []
# create a list of counts for non-cnv genes parallel to labelnames list
non_cnv_genes = []

# loop over species in name
for sp in names:
    cnv_genes.append(gene_counts[sp][1])
    non_cnv_genes.append(gene_counts[sp][2])
    
# create figure
fig = plt.figure(1, figsize = (4.3,2.56))

# add axe to fig
ax = fig.add_subplot(1, 1, 1)

# Set the bar width
bar_width = 0.75

# set positions of the left bar-boundaries
bar_left = [i+1 for i in range(len(cnv_genes))]

# set positions of the x-axis ticks (center of the bars as bar labels)
tick_pos = [i+(bar_width/2) for i in bar_left]

# Create a bar plot, in position bar_left for cnv genes
ax.bar(bar_left, cnv_genes, width=bar_width, label = 'CNV', color= 'black')

# Create a bar plot, in position bar_left for non_cnv genes on top of cnv_genes
ax.bar(bar_left, non_cnv_genes, width=bar_width, bottom= cnv_genes, label = 'non-CNV', color = 'white')

# set the x ticks with names
if species == 'species':
    plt.xticks(tick_pos, labelnames, style = 'italic')
elif species == 'human':
    plt.xticks(tick_pos, labelnames)

# set axis labels
ax.set_ylabel('Number of genes\nwith short 3\'UTR', size = 12, ha = 'center', fontname = 'Arial', family = 'sans-serif', color = 'grey')

# check if species is species or human
if species == 'human':
    # add a x-axis label
    ax.set_xlabel('DGV releases', size = 12, ha = 'center', fontname = 'Arial', family = 'sans-serif')
elif species == 'species':
    # add a x-axis label
    ax.set_xlabel('Vertebrate species', size = 12, ha = 'center', fontname = 'Arial', family = 'sans-serif')


# Set a buffer around the edge
plt.xlim([min(tick_pos)-bar_width, max(tick_pos)+bar_width])



# add a light grey horizontal grid to the plot, semi-transparent, 
ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)  
# hide these grids behind plot objects
ax.set_axisbelow(True)

# do not show lines around figure  
ax.spines["top"].set_visible(False)    
ax.spines["bottom"].set_visible(False)    
ax.spines["right"].set_visible(False)    
ax.spines["left"].set_visible(False)  


# do not show ticks
  
plt.tick_params(
    axis='y',       # changes apply to the x-axis and y-axis (other option : x, y)
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    right = 'off',
    left = 'off',          
    labelbottom='off', # labels along the bottom edge are off 
    colors = 'grey')  



#
# remove top axes and right axes ticks
ax.get_xaxis().tick_bottom()
#ax.get_yaxis().tick_left()

# get outputfile
outputfile = 'Fig_short3UTR_counts_' + '_' + species + '_' + cnv_length + '_' + chromos
print(outputfile)
  
# save figure
fig.savefig(outputfile + '.eps', bbox_inches = 'tight')
    