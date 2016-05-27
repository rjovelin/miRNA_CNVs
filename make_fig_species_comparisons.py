# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 10:21:42 2016

@author: RJovelin
"""

# make figure comparing normalized sites in 5'UTR or CDS between human and each other 
# species for genes that have short 3'UTRs 

import os
import sys
import matplotlib.pyplot as plt
import numpy as np

# usage python3 make_fig_species_comparisons.py inputfile

# get filename from command
filename = sys.argv[1]

# get domain from filename
if '5UTR' in filename:
    domain  = '5UTR'
elif 'CDS' in filename:
    domain = 'CDS'
    
# get the gene types (CNV or not_CNV) from filename
if 'not_CNV' in filename:
    gene_type = 'not_CNV'
else:
    gene_type = 'CNV'
    
# check if all chromos (including unplaced, unlocated, and MT) are used
# or if only valid chromos are used 
# note: only CNVs on valid chromos are reported in DGV, so if all chromos are
# used it may introduce a bias by calling non CNV genes genes that cannot be detected

if 'valid_chromos' in filename:
    chromos = 'valid_chromos'
elif 'all_chromos' in filename:
    chromos = 'all_chromos'

# get the type of CNVs to consider from the filename
if 'all_CNVs' in filename:
    cnv_length = 'CNV_all_length'
elif 'long_CNVs' in filename:
    cnv_length = 'CNV_greater_1Kb'

# get DGV release version
release_version = filename[filename.index('GRCh37') : filename.rfind('_' + gene_type)]

print(domain, gene_type, chromos, cnv_length, release_version)


# create a dict to store data {Hsa_sp2: [[Hsa_sites], [sp2_sites]]} 
species_data = {}

# parse inputfile
infile = open(filename, 'r')
# skip header
infile.readline()

# loop over file
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split('\t')
        # get comp_names
        comp_names = line[0]
        # get species 
        species = line[1]
        # check if species pairs are already in dict
        if comp_names not in species_data:
            # double check that species is human
            assert species == 'Hsa', 'species is not human'
            # populate dict, initialize empty list
            species_data[comp_names] = []
        elif comp_names in species_data:
            # double check that species is not human
            assert species != 'Hsa', 'species is human'
        # convert values to float, add values to list
        sites = list(map(lambda x: float(x), line[2:]))
        # copy sites to list
        species_data[comp_names].append(sites[:])

# close file after reading
infile.close()

# loop over comp_names, print comp_names and number of sites for each species
for comp_names in species_data:
    print(comp_names, len(species_data[comp_names][0]), len(species_data[comp_names][1]))
    print(comp_names, np.mean(species_data[comp_names][0]), np.mean(species_data[comp_names][1]))

# make a list of species pairs
species_pairs = ['Hsa_Ptr', 'Hsa_Mml', 'Hsa_Mmu', 'Hsa_Bta', 'Hsa_Gga']


# make a list of lists with all data
all_data = []
# loop over species pairs, append list of sites
for i in range(len(species_pairs)):
    all_data.append(species_data[species_pairs[i]][0])
    all_data.append(species_data[species_pairs[i]][1])
print('data consolidated in array')
    
    
# create figure
fig = plt.figure(1, figsize = (4.3,2.56))

# add a plot to figure (1 row, 1 column, 1 plot)
ax = fig.add_subplot(1, 1, 1)    

# write label for y axis
ytext = ax.set_ylabel('Normalized number of miRNA\nsites per gene', color = 'grey', size = 10, ha = 'center', fontname = 'Arial', family = 'sans-serif')

# set tick label
names = []
for i in range(len(species_pairs)):
    names.extend(list(species_pairs[i].split('_')))


# add labels to x-ticks, and align center, set size to 10
ax.set_xticklabels(names, ha = 'center', size = 10, fontname = 'Arial', family = 'sans-serif')


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
    
## note some ways to modify label font size and color of the ticks
#for label in ax.get_yticklabels():
#    label.set_color('black')
#    label.set_fontsize(10)

# add title
if gene_type == 'CNV':
    plt.title('CNV genes', size = 10, fontname = 'Arial')  
elif gene_type == 'not_CNV':
    plt.title('non-CNV genes', size = 10, fontname = 'Arial')

# get outputfile
outputfile = 'Fig_species_comp_' + domain + '_' + chromos + '_' + cnv_length + '_' + release_version + '_' + gene_type
print(outputfile)

  
# save figure
fig.savefig(outputfile + '.eps', bbox_inches = 'tight')
    
    
    
    
    
    
   