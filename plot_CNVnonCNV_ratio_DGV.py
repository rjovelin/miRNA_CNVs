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

  