# -*- coding: utf-8 -*-
"""
Created on Fri Jul 22 17:40:04 2016

@author: RJovelin
"""

# use this script to plot the number of short 3'UTR genes in CNV and non-CNVs for release of the DGV

# usage PlotNumberShort3UTRGenes.py [options]
# [7/15]: 3'UTR length, genes with 3'UTR length < 7bp or < 15bp are considered short 3'UTR genes

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


# use all chromos (including unplaced, unlocated, and MT) or only valid chromos 
# note: only CNVs on valid chromos are reported in DGV, so if all chromos are
# used it may introduce a bias by calling non CNV genes genes that cannot be detected

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

# get the minimum 3'UTR length
L = int(sys.argv[1])
assert L in [15, 7], 'minimum 3UTR length is not correct'
print(L)

# Count the number of short 3'UTR genes in CNV and non-CNV in each release
# create a dict {release: [N short 3'UTR, N short 3UTR CNV genes, N short 3UTR non-CNV genes]}
ShortGenes = {}

# get UTR file
UTR_file = 'H_sapiens_3UTR_length_' + chromos + '.txt'
print(UTR_file)

# make a list of CNV files
CNV_files = ['H_sapiens_GRCh37_2013-05_' + cnv_length + '_' + chromos + '.txt',
             'H_sapiens_GRCh37_2013-07_' + cnv_length + '_' + chromos + '.txt',
             'H_sapiens_GRCh37_2014_' + cnv_length + '_' + chromos + '.txt',
             'H_sapiens_GRCh37_2015_' + cnv_length + '_' + chromos + '.txt']

# loop over CNV files
for filename in CNV_files:
    print(filename)
    # sort genes based on CNV status
    CNV_status = sort_genes_CNV_status(filename)
    print(len(CNV_status))
    # sort genes based on 3' UTR length
    UTR_length = sort_genes_3UTR_length(UTR_file, L)
    print(len(UTR_length))
    # get release version
    release_version = filename[filename.index('GRCh37'): filename.index('_CNV')]
    print(release_version)
    # count total number of genes with short 3'UTR
    total_short = 0
    # count CNV genes with short UTR
    cnv_short = 0
    # count non-CNV genes with short UTR
    non_cnv_short = 0    
    # loop over genes in UTR_length
    for gene in UTR_length:
        if UTR_length[gene] == 'short':
            total_short += 1
            # check CNV status
            if CNV_status[gene] == 'CNV':
                cnv_short += 1
            elif CNV_status[gene] == 'not_CNV':
                non_cnv_short += 1
            
    # check that numbers add up
    assert total_short == cnv_short + non_cnv_short, 'sum cnv and non-cnv short is not equal to total short'
    # populare dict
    ShortGenes[release_version] = [total_short, cnv_short, non_cnv_short]
print('got short gene counts for each species')


# plot the number of genes with short 3'UTR in each release of the DVG as a bar graph

Releases = ['GRCh37_2013-05', 'GRCh37_2013-07', 'GRCh37_2014', 'GRCh37_2015']
labelnames = ['2013a', '2013b', '2014', '2015']
    
# create list of counts for cnv genes and non-CNV genes parallel to labelnames list
cnv_genes, non_cnv_genes = [], []
for i in Releases:
    cnv_genes.append(ShortGenes[i][1])
    non_cnv_genes.append(ShortGenes[i][2])


###########################


# count short 3'UTR genes for each study of each release of the DGV

# create a dict {release: {study: [pubmedID, total, CNV, non-CNV]}}
StudiesShortGenes = {}


# loop over CNV files
for filename in CNV_files:
    print(filename)
    # sort genes based on CNV status
    CNV_status = sort_genes_CNV_status(filename)
    print(len(CNV_status))
    # sort genes based on 3' UTR length
    UTR_length = sort_genes_3UTR_length(UTR_file, L)
    print(len(UTR_length))
    # get release version
    release_version = filename[filename.index('GRCh37'): filename.index('_CNV')]
    print(release_version)
    
    
####################### continue here









    



# write headers to file
newfile.write('# Number of Human genes with short (< 7 bp) 3\'UTR\n')
newfile.write('\t'.join(['Study', 'PubMedID', 'Total', 'CNV', 'non-CNV']) + '\n')

# get the dictionaries of {reference: pubmedid} for studies reported in the CGV
references = get_DGV_references(CNV_file)

# get synonym names for all genes {gene name : [list of synonyms]}
synonyms = get_synonyms('H_sapiens.gff3')
    
# loop over study
for study in references:
    print(study)
    # get the set of CNV genes corresponding to that study
    CNV_genes = get_human_CNV_genes_single_study(CNV_file, study, CNV_size)
    
    print('# CNV genes', len(CNV_genes))    
    
    # sort genes based on 3' UTR length
    UTR_length = sort_genes_3UTR_length(UTR_file)
    print(len(UTR_length))

    # get the CNV status of all short #' UTR genes
    # create a dict {gene: CNV status}    
    CNV_status = {}
    
    # loop over gene in UTR_length
    for gene in UTR_length:
        # set boolean
        is_cnv = False
        # check UTR length
        if UTR_length[gene] == 'short':
            # ask if gene in CNV genes
            if gene in CNV_genes or gene.upper() in CNV_genes:
                # gene is CNV, add gene and status to dict
                CNV_status[gene] = 'CNV'
            else:
                # ask if any of the gene synonyms are in CNV genes
                for name in synonyms[gene]:
                    # check if in CNV genes
                    if name in CNV_genes or name.upper() in CNV_genes:
                        # update boolean variable
                        is_cnv = True
                # check if gene is CNV
                if is_cnv == True:
                    CNV_status[gene] = 'CNV'
                elif is_cnv == False:
                    CNV_status[gene] = 'not_CNV'
    
    print('# short genes with CNV status', len(CNV_status))
    
    # count total number of genes with short 3'UTR
    total_short = 0
    # count CNV genes with short UTR
    cnv_short = 0
    # count non-CNV genes with short UTR
    non_cnv_short = 0    
    
    # loop over genes in UTR_length
    for gene in UTR_length:
        if UTR_length[gene] == 'short':
            total_short += 1
            # check CNV status
            if CNV_status[gene] == 'CNV':
                cnv_short += 1
            elif CNV_status[gene] == 'not_CNV':
                non_cnv_short += 1
            
    print('total short', total_short)
    print('CNV short', cnv_short)
    print('non_CNV', non_cnv_short)
    
    assert total_short == cnv_short + non_cnv_short, 'sum cnv and non-cnv short is not equal to total short'
        
    # write results to file
    newfile.write('\t'.join([study, references[study], str(total_short), str(cnv_short), str(non_cnv_short)]) + '\n')
    
    
# close file after writing
newfile.close()

Contact GitHub API Training Shop Blog About
© 2016 GitHub, Inc. Terms Privacy Security Status Help


#######################



# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 23:45:20 2015

@author: Richard
"""




# plot the number of studies with ratio of CNV / non-CNV short 3'UTR genes
# for each release of the DGV

# usage plot_CNVnonCNV_ratio_DGV.py [parameters]
# [True/False] use valid chromos or all chromos
# [all_CNVs/long_CNVs] use all CNVs or CNVs > 1Kb


import os
import sys
import matplotlib.pyplot as plt

# get chromos from command
keep_valid_chromos = sys.argv[1]
if keep_valid_chromos == 'True':
    chromos = 'valid_chromos'
elif keep_valid_chromos == 'False':
    chromos = 'all_chromos'
    
# get the type of CNVs to consider from the command [long_CNVs or all_CNVs]
CNV = sys.argv[2]
if CNV == 'all_CNVs':
    cnv_length = 'CNV_all_length'
elif CNV == 'long_CNVs':
    cnv_length = 'CNV_greater_1Kb'


# make a list of files with number of genes with short 3' UTR
# in CNV and non-CNV for each version of DGV

files = [i for i in os.listdir() if ('Human_Counts_Short3UTR_' + cnv_length + '_' + chromos) in i]

#create a dict with version and list of counts
ratio = {}
for filename in files:
    # grab version
    version = filename[filename.rindex('_')+1: filename.index('.txt')]
    # open file for reading
    infile = open(filename, 'r')
    # skip 2 first lines
    infile.readline()
    infile.readline()
    # create a list of 0
    counts = [0] * 10
    # loop over file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split('\t')
            cnv = int(line[3])
            noncnv = int(line[4])
            freq = cnv / noncnv
            # get the index in list where value should go
            pos = int((freq * 100) // 10)
            counts[pos] += 1
    # close file
    infile.close()
    # populate dict with version : empty list pairs
    ratio[version] = counts
    
# ~1.33x wider than tall
# Common sizes: (10, 7.5) and (12, 9)    
plt.figure(figsize=(8, 5.5))    
  
# Remove the plot frame lines.     
ax = plt.subplot(111)    
ax.spines["top"].set_visible(False)    
ax.spines["bottom"].set_visible(True)    
ax.spines["right"].set_visible(False)    
ax.spines["left"].set_visible(True)    

# find the maximum y value
maximum = 0
for year in ratio:
    for val in ratio[year]:
        if val > maximum:
            maximum = val

# Limit the range of the plot to data
plt.ylim(0, maximum)    
plt.xlim(0, 10)    
  
# adjust size of ticks    
  # set major ticks on the y axis
plt.yticks(range(0, maximum + 10, 10), fontname = 'Arial', fontsize=16)    
plt.xticks([i for i in range(0, 11)], fontname = 'Arial', fontsize=16)    
  
# add a light grey horizontal grid to the plot, semi-transparent, 
ax.yaxis.grid(True, linestyle='-', which='major', color='0.80', alpha=0.6)  
# hide these grids behind plot objects
ax.set_axisbelow(True)

# Remove the tick marks; they are unnecessary with the tick lines we just plotted.    
plt.tick_params(axis="both", which="both", bottom="on", top="off",    
                labelbottom="on", left="on", right="off", labelleft="on")    

# plot the data

releases = ['2013a', '2013b', '2014', '2015']

# make a list of years
years = ['2013-05', '2013-07', '2014', '2015']

for year in years:
    if year == '2013-05':
        plt.plot([i + 0.5 for i in range(0, 10)], ratio[year], linestyle = '-', color = '0.80', marker = 'o', markersize = 10, markeredgewidth = 2, markerfacecolor = '0.80', markeredgecolor = '0.80', lw = 3, label = '2013a')
    elif year == '2013-07':
        plt.plot([i + 0.5 for i in range(0, 10)], ratio[year], linestyle = '-', color = '0.60', marker = 'o', markersize = 10, markeredgewidth = 2, markerfacecolor = '0.60', markeredgecolor = '0.60',   lw = 3, label = '2013b')
    elif year == '2014':
        plt.plot([i + 0.5 for i in range(0, 10)], ratio[year], linestyle = '-', color = '0.40', marker = 'o', markersize = 10, markeredgewidth = 2, markerfacecolor = '0.40', markeredgecolor = '0.40', lw = 3, label = '2014')
    elif year == '2015':
        plt.plot([i + 0.5 for i in range(0, 10)], ratio[year], linestyle = '-', color = '0.20', marker = 'o', markersize = 10, markeredgewidth = 2, markerfacecolor = '0.20', markeredgecolor = '0.20', lw = 3, label = '2015')

# add x and y labels
plt.xlabel('Ratio of number of CNV genes / non-CNV genes', fontname = 'Arial', fontsize = 16)
plt.ylabel('Number of studies in DGV', fontname = 'Arial', fontsize = 16)


# do not show ticks
plt.tick_params(
    axis='y',       # changes apply to the x-axis and y-axis (other option : x, y)
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    right = 'off',
    left = 'off',          
    labelbottom='off', # labels along the bottom edge are off 
    colors = 'grey'
    )  

# do not show lines around figure  
ax.spines["top"].set_visible(False)    
ax.spines["bottom"].set_visible(False)    
ax.spines["right"].set_visible(False)    
ax.spines["left"].set_visible(False) 

# get maximum y value
ymax = 0
for year in ratio:
    if max(ratio[year]) > ymax:
        ymax = max(ratio[year])
        print(ymax)

# Set a buffer around the edge
plt.ylim([-1, round(ymax, -1)])


# bbox_inches="tight" removes all the extra whitespace on the edges of the plot    
plt.savefig("Fig_ratio_CNVnonCNV_DGV.eps", bbox_inches="tight")  

  























######################################


# create figure
fig = plt.figure(1, figsize = (4.3,2.56))

# add axe to fig
ax = fig.add_subplot(1, 1, 1)

# Set the bar width
bar_width = 0.5

# set positions of the x-axis ticks
xtickpos = [0, 0.7, 1.4, 2.1]

# Create a bar plot for cnv genes
ax.bar(xtickpos, cnv_genes, width=bar_width, label = 'CNV', color= '#ef8a62')
# Create a bar plot for non_cnv genes on top of cnv_genes
ax.bar(xtickpos, non_cnv_genes, width=bar_width, bottom= cnv_genes, label = 'non-CNV', color = '#67a9cf')

# set font for all text in figure
FigFont = {'fontname':'Arial'} 

# add labels to x axis ticks
plt.xticks(xtickpos, labelnames, **FigFont)

# set axis labels
ax.set_ylabel('Number of genes\nwith short 3\'UTR', size = 10, ha = 'center', color = 'black', **FigFont)

# do not show lines around figure  
ax.spines["top"].set_visible(False)    
ax.spines["bottom"].set_visible(True)    
ax.spines["right"].set_visible(False)    
ax.spines["left"].set_visible(False)  
# offset the spines
for spine in ax.spines.values():
    spine.set_position(('outward', 5))

# add a light grey horizontal grid to the plot, semi-transparent, 
ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)  
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
    labelsize = 10,
    direction = 'out') # ticks are outside the frame when bottom = 'on'  
      
# Set the tick labels font name
for label in ax.get_yticklabels():
    label.set_fontname('Arial')
    
# create a margin around the x axis
plt.margins(0.05)

# get outputfile
if L == 7:
    outputfile = 'PlotShort3UTRCountsDGV_7bp_' + cnv_length + '_' + chromos
elif L == 15:
    outputfile = 'PlotShort3UTRCountsDGV_15bp_' + cnv_length + '_' + chromos
print(outputfile)
  
# save figure
fig.savefig(outputfile + '.eps', bbox_inches = 'tight')
