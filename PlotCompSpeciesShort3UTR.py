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
    
        
################################################
    

# plot the number of target sites of mirna families shared between human and
# with each other species for CNV and for non-CNV genes
 

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
outputfile = 'truc_' + domain + '_' + chromos + '_' + cnv_length + '_' + release_version + '_' + gene_type
print(outputfile)

  
# save figure
fig.savefig('truc.pdf', bbox_inches = 'tight')



