# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 12:12:16 2015

@author: RJovelin
"""

# plot the frequency of studies for which the number of mirna sites is greater
# for CNV genes, lower for CNV genes and for which there is no significant 
# difference between CNV and non-CNv genes


# usage python3 plot_freq_individual_studies.py options
# [targetscan/miranda] predictor used to predict mirna target sites
# outputfile

import os
import sys
import matplotlib.pyplot as plt
from matplotlib import patches as mpatches
import numpy as np


# get the predictor [targetscan or miranda] 
predictor = sys.argv[1]
print(predictor)

# get outputfile
outputfile = sys.argv[2]
print(outputfile)

# create a list of summary files for each release of the DGV 
files = [i for i in os.listdir() if 'H_sapiens_single' in i and predictor in i]
# sort filenames
files.sort()

# create a dict {DGV_release : [N_cnv_greater, N_cnv_lower, N_no_diff]}
studies = {}

# loop over filename
for filename in files:
    print(filename)
    # open file for reading
    infile = open(filename, 'r')
    # skip header
    header = infile.readline()
    # count number of studies with CNV greater, CNV lower and no-diff
    CNV_greater, CNV_lower, no_diff = 0, 0, 0
    # loop over file
    for line in infile:
        # consider only normalized binding sites
        line = line.rstrip()
        if line != '':
            line = line.split('\t')
            # get CNV normalized targets
            CNV_targets = float(line[8])
            # get non-CNV normalized targets
            nonCNV_targets = float(line[10])
            # get P value
            P_val = float(line[12])
            # check if difference is significant
            if P_val > 0.05:
                # no significant difference
                no_diff += 1
            elif P_val < 0.05:
                # check if CNV greater or lower
                if CNV_targets > nonCNV_targets:
                    CNV_greater += 1
                elif CNV_targets < nonCNV_targets:
                    CNV_lower += 1
    # close file
    infile.close()
    # populate dict with frequencies
    # get total number of studies for given release
    total = CNV_greater + CNV_lower + no_diff
    
    if '2013-05' in filename:
        studies['2013a'] = [CNV_greater / total, CNV_lower / total, no_diff / total]
    elif '2013-07' in filename:
        studies['2013b'] = [CNV_greater / total, CNV_lower / total, no_diff / total]
    elif '2014' in filename:
        studies['2014'] = [CNV_greater / total, CNV_lower / total, no_diff / total]
    elif '2015' in filename:
        studies['2015'] = [CNV_greater / total, CNV_lower / total, no_diff / total]

for i in studies:
    print(i, studies[i])        

# make a list of release names
releases = ['2013a', '2013b', '2014', '2015']

# create parallel lists with CNV_greater, CNV_lower and no_diff for each release
cnv_greater, cnv_lower, nodiff = [], [], []

# loop over each release
for i in releases:
    print(i, sum(studies[i]))
    # get the corresponding frequencies
    cnv_greater.append(studies[i][0])
    cnv_lower.append(studies[i][1])
    nodiff.append(studies[i][2])

# make a list with added values between cn_greater and cnv_lower
cnv_added = []
for i in range(len(cnv_greater)):
    cnv_added.append(cnv_greater[i] + cnv_lower[i])

  
# create figure
fig = plt.figure(1, figsize = (4.3,2.56))

# add axe to fig
ax = fig.add_subplot(1, 1, 1)

# Set the bar width
bar_width = 0.75

# set positions of the left bar-boundaries
bar_left = [i+1 for i in range(len(cnv_greater))]

# set positions of the x-axis ticks (center of the bars as bar labels)
tick_pos = [i+(bar_width/2) for i in bar_left]

# Create a bar plot, in position bar_left for cnv greater
plt.bar(bar_left, cnv_greater, width=bar_width, color= 'black')

# Create a bar plot, in position bar_left for cnv lower on top of cnv_greater
plt.bar(bar_left, cnv_lower, width=bar_width, bottom= cnv_greater, color = 'grey')

# create a bar plot, in position bar_left for no diff on top of cnv lower
plt.bar(bar_left, nodiff, width = bar_width, bottom = cnv_added , color = 'white')

# set the x ticks with names
plt.xticks(tick_pos, releases, size = 12)

# set the y ticks
plt.yticks([i/100 for i in range(0, 125, 25)], [0, 0.25, 0.50, 0.75, 1])

# set axis labels
plt.ylabel('Frequency of studies in DGV', size = 12, ha = 'center', fontname = 'Helvetica', family = 'sans-serif', color = 'grey')

plt.xlabel('DGV releases', size = 12, ha = 'center', fontname = 'Helvetica', family = 'sans-serif')

# Set a buffer around the edge
plt.xlim([min(tick_pos)-bar_width, max(tick_pos)+bar_width])

# add a light grey horizontal grid to the plot, semi-transparent, 
ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)  
# hide these grids behind plot objects
ax.set_axisbelow(True)


plt.margins()
  
# do not show lines around figure  
ax.spines["top"].set_visible(False)    
ax.spines["bottom"].set_visible(False)    
ax.spines["right"].set_visible(False)    
ax.spines["left"].set_visible(False)      

if predictor == 'targetscan':
    plt.title('TargetScan', size = 12)  
elif predictor == 'miranda':
    plt.title('miRanda', size = 12)
  
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







# remove top axes and right axes ticks
#ax.get_xaxis().tick_bottom()
#ax.get_yaxis().tick_left()


  
# save figure
fig.savefig(outputfile, bbox_inches = 'tight')
    