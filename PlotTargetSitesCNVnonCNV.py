# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 15:16:57 2016

@author: RJovelin
"""

# use this script to compare target sites and normalized target sites between CNV and non-CNv genes in each species



# usage python3 PlotTargetSitesCNVnonCNV.py options
# [3UTR/5UTR/CDS]: region to consider
# [True/False]
# [long_CNVs/all_CNVs]





from CNV_miRNAs import *
import os
import sys
import matplotlib.pyplot as plt










# get the region to consider to predict target sites [3UTR or 5UTr or CDS]
domain = sys.argv[1]
print(domain)

 

# get the option to keep genes on all chromos (False) or only on assembled nuclear chromosomes (True) 
keep_valid_chromos = sys.argv[3]
if keep_valid_chromos == 'True':
    chromos = 'valid_chromos'
elif keep_valid_chromos == 'False':
    chromos = 'all_chromos'
print(keep_valid_chromos, chromos)

# get the type of CNVs to consider from the command [long_CNVs or all_CNVs]
long_CNV = sys.argv[4]
if long_CNV == 'all_CNVs':
    cnv_length = 'CNV_all_length'
elif long_CNV == 'long_CNVs':
    cnv_length = 'CNV_greater_1Kb'
print(long_CNV, cnv_length)

# get the type of data to consider [targets or normalized or length]
data_type = sys.argv[5]

# make a list of species names
species_names = ['H_sapiens', 'P_troglodytes', 'M_mulatta', 'M_musculus', 'B_taurus', 'G_gallus']

species_codes = {'H_sapiens': 'Hsa', 'P_troglodytes': 'Ptr', 'M_mulatta': 'Mml',
                 'M_musculus': 'Mmu', 'B_taurus': 'Bta', 'G_gallus':'Gga'}

# create a dict of {species name : [[CNV], [non-CNV]]}
species_data = {}

# loop over species names
for species in species_names:
    print(species)
    
    # get the summary table
    if species == 'H_sapiens':
        # get summary table for DGV 2015 release
        summary_table = 'H_sapiens_' + domain + '_summary_' + predictor + '_' + chromos + '_' + cnv_length + '_GRCh37_2015.txt'
    else:
        # get summary table 
        summary_table = species + '_' + domain + '_summary_' + predictor + '_' + chromos + '_' + cnv_length + '.txt'
    
    
    # parse the summary table into a list of list of values
    # [CNV_targets, nonCNV_targets, CNV_seq, nonCNV_seq, CNV_normalized, nonCNV_normalized]    
    parsed_data = parse_summary_table_targets(summary_table)
    
    # populate dicts for the type of data to consider
    if data_type == 'targets':
        species_data[species] = [parsed_data[0], parsed_data[1]]
    elif data_type == 'normalized':
        species_data[species] = [parsed_data[4], parsed_data[5]]
    elif data_type == 'length':
        species_data[species] = [parsed_data[2], parsed_data[3]]        
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
# check data type
if data_type == 'targets':
    ytext = ax.set_ylabel('Number of miRNA sites per gene', color = 'grey', size = 10, ha = 'center', fontname = 'Arial', family = 'sans-serif')
elif data_type == 'normalized':
    ytext = ax.set_ylabel('Normalized number of miRNA\nsites per gene', color = 'grey', size = 10, ha = 'center', fontname = 'Arial', family = 'sans-serif')
elif data_type == 'length':
    ytext = ax.set_ylabel('Sequence length', color = 'grey', size = 10, ha = 'center', fontname = 'Arial', family = 'sans-serif')

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
    
## note some ways to modify label font size and color of the ticks
#for label in ax.get_yticklabels():
#    label.set_color('black')
#    label.set_fontsize(10)

# add title
if predictor == 'targetscan':
    plt.title('TargetScan', size = 10, fontname = 'Arial')  
elif predictor == 'miranda':
    plt.title('miRanda', size = 10, fontname = 'Arial')

# get outputfile
outputfile = 'Fig_sites_CNVvsNonCNV_' + domain + '_' + chromos + '_' + cnv_length + '_' + predictor + '_' + data_type
print(outputfile)

  
# save figure
fig.savefig(outputfile + '.eps', bbox_inches = 'tight')
    
    
    
    
    
    
   