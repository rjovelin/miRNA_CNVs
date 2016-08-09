# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 11:59:52 2016

@author: RJovelin
"""


# use this script to compare the number of miRNA targets per nucleotide
# between human and each other vertebrate species for genes with short 3'UTRs


# usage PlotCompSpeciesShort3UTR.py [options]
# - [targetscan/miranda]: algorithm used to predict target sites
# - [5UTR/CDS]: domain to compare


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




# get the prediction method
predictor = sys.argv[1]
# get the region to consider to predict target sites [3UTR or 5UTr or CDS]
domain = sys.argv[2]
assert domain in ['5UTR', 'CDS'], 'the gene domain is should be 5UTR or CDS'

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



















from CNV_miRNAs import *
import os
import sys
import numpy as np
from scipy import stats
import math


# usage python3 interspecies_comparison_targets_short3UTR_genes.py options
# [5UTR/CDS]
# [targetscan/miranda]
# [True/False]
# [long_CNVs/all_CNVs]
# [shared_miRs/all_miRs]
# human_CNV_file 

# get the region to consider to predict target sites [5UTr or CDS]
domain = sys.argv[1]
print(domain)

# get the predcitor used to predictor target sites [targetscan or miranda]
predictor = sys.argv[2]
print(predictor)

# get the option to keep genes on all chromos (False) or only on assembled 
# nuclear chromosomes only from the command
keep_valid_chromos = sys.argv[3]

# check if all chromos (including unplaced, unlocated, and MT) are used
# or if only valid chromos are used 
# note: only CNVs on valid chromos are reported in DGV, so if all chromos are
# used it may introduce a bias by calling non CNV genes genes that cannot be detected

if keep_valid_chromos == 'True':
    keep_valid_chromos = True
    chromos = 'valid_chromos'
elif keep_valid_chromos == 'False':
    keep_valid_chromos = False
    chromos = 'all_chromos'
print(keep_valid_chromos, chromos)


# get the type of CNVs to consider from the command
long_CNV = sys.argv[4]
if long_CNV == 'all_CNVs':
    cnv_length = 'CNV_all_length'
elif long_CNV == 'long_CNVs':
    cnv_length = 'CNV_greater_1Kb'
print(long_CNV, cnv_length)

# get the option to consider all miRNAs or conserved miRNAs between species pairs
# [shared_miRs or all_miRs]
conservation = sys.argv[5]
print(conservation)

# get the human CNV file from the command, so that comparisons can be made 
# between all species and different version of the human DVG
human_CNV_file = sys.argv[6]
print(human_CNV_file)

# get realease version of the DGV
release_version = human_CNV_file[human_CNV_file.index('GRCh37') : human_CNV_file.index('_CNV')]

# get outputfile
outputfile = 'Interspecies_comp_targets_' + predictor + '_' + domain + '_short3UTR_' + chromos + '_' + cnv_length + '_' + conservation + '_' + release_version + '.txt'
print(outputfile)

# open file for writing
newfile = open(outputfile, 'w')

# write header to file
newfile.write('Species1' + '\t' + 'Species2' + '\t')

# write number of miRNAs if using shared miRNAs
if conservation  == 'shared_miRs':
    newfile.write('shared_miRNAs' + '\t')

newfile.write('\t'.join(['N_CNV_genes_sp1', 'N_CNV_genes_sp2',
                         'CNV_mean_targets_sp1', 'CNV_SEM_targets_sp1',
                         'CNV_mean_targets_sp2', 'CNV_SEM_targets_sp2', 'P_diff_CNV_targets',
                         'N_nonCNV_genes_sp1', 'N_nonCNV_genes_sp2', 
                         'nonCNV_mean_targets_sp1', 'nonCNV_SEM_targets_sp1',
                         'nonCNV_mean_targets_sp2', 'nonCNV_SEM_targets_sp2', 'P_diff_nonCNV_targets',
                         'CNV_mean_normalized_targets_sp1', 'CNV_SEM_normalized_targets_sp1',
                         'CNV_mean_normalized_targets_sp2', 'CNV_SEM_normalized_targets_sp2', 'P_diff_norm_CNV_targets',
                         'nonCNV_mean_normalized_targets_sp1', 'nonCNV_SEM_normalized_targets_sp1',
                         'nonCNV_mean_normalized_targets_sp2', 'nonCNV_SEM_normalized_targets_sp2', 'P_diff_norm_nonCNV_targets',
                         'CNV_mean_seq_length_sp1', 'CNV_SEM_seq_length_sp1',
                         'CNV_mean_seq_length_sp2', 'CNV_SEM_seq_length_sp2', 'P_diff_seq_CNV_length',
                         'nonCNV_mean_seq_length_sp1', 'nonCNV_SEM_seq_length_sp1',
                         'nonCNV_mean_seq_length_sp2', 'nonCNV_SEM_seq_length_sp2', 'P_diff_seq_nonCNV_length']) + '\n')                        

# create a lamda function to transform value into string
Gstr = lambda x: str(x)

# make a list of species names
species = ['H_sapiens', 'P_troglodytes', 'M_mulatta', 'M_musculus', 'R_norvegicus', 'B_taurus', 'C_familiaris', 'G_gallus']


# loop over species
for i in range(len(species) -1):
    # compare pairs of species
    for j in range(i+1, len(species)):
        print(species[i], species[j])
        # get genome files        
        genome_sp1 = species[i] + '_genome.txt'
        genome_sp2 = species[j] + '_genome.txt'
        # get chromos files
        valid_chromos_sp1 = species[i] + '_valid_chromos.txt'
        valid_chromos_sp2 = species[j] + '_valid_chromos.txt'
        # get GFF annotation files
        GFF_sp1 = species[i] + '.gff3'
        GFF_sp2 = species[j] + '.gff3'
        print(genome_sp1, genome_sp2)
        print(valid_chromos_sp1, valid_chromos_sp2)
        print(GFF_sp1, GFF_sp2)
        
        # get targetscan sequence input file
        seq_input_sp1 = species[i] + '_' + domain + '_' + chromos + '_targetscan.txt'
        seq_input_sp2 = species[j] + '_' + domain + '_' + chromos + '_targetscan.txt'
        print(seq_input_sp1, seq_input_sp2)
        
        # get predictor output
        predicted_targets_sp1 = species[i] + '_' + domain + '_' + chromos + '_predicted_sites_' + predictor + '.txt'
        predicted_targets_sp2 = species[j] + '_' + domain + '_' + chromos + '_predicted_sites_' + predictor + '.txt'
        print(predicted_targets_sp1, predicted_targets_sp2)        
                
        # get UTR files
        UTR_sp1 = species[i] + '_3UTR_length_' + chromos + '.txt'      
        UTR_sp2 = species[j] + '_3UTR_length_' + chromos + '.txt'
        print(UTR_sp1, UTR_sp2)
         
        # get CNV files
        if species[i] == 'H_sapiens':
            CNV_file_sp1 = human_CNV_file
            CNV_file_sp2 = species[j] + '_' + cnv_length + '_' + chromos + '.txt'
        else:
            CNV_file_sp1 = species[i] + '_' + cnv_length + '_' + chromos + '.txt'
            CNV_file_sp2 = species[j] + '_' + cnv_length + '_' + chromos + '.txt'
        print(CNV_file_sp1, CNV_file_sp2)
        
        # get mature mirna files
        mature_sp1 = species[i] + '_mature.txt'
        mature_sp2 = species[j] + '_mature.txt'        
        print(mature_sp1, mature_sp2)
        
        # get a set of shared seeds between species pair
        shared_seeds = find_conserved_mirna_families(mature_sp1, mature_sp2)
        
        # get CNV gene status
        CNV_status1 = sort_genes_CNV_status(CNV_file_sp1)
        CNV_status2 = sort_genes_CNV_status(CNV_file_sp2)
        print('CNV status', len(CNV_status1), len(CNV_status2))
        
        # parse predictor outputfile {gene: [N_sites, seq_length, N_sites/seq_length]}
        # check predictor
        if predictor == 'targetscan':
            # check if all mirna families are considered or only families shared by the 2 species
            if conservation == 'shared_miRs':
                # filter targets for shared mirnas
                targets_sp1 = parse_targetscan_output(seq_input_sp1, predicted_targets_sp1, 'conserved', shared_seeds, mature_sp1)
                targets_sp2 = parse_targetscan_output(seq_input_sp2, predicted_targets_sp2, 'conserved', shared_seeds, mature_sp2)
            elif conservation == 'all_miRs':
                # take all sites for all mirnas
                targets_sp1 = parse_targetscan_output(seq_input_sp1, predicted_targets_sp1, 'all')
                targets_sp2 = parse_targetscan_output(seq_input_sp2, predicted_targets_sp2, 'all')
        
        elif predictor == 'miranda':
            # check if all mirna families are considered or only families shared by the 2 species
            if conservation == 'shared_miRs':
                targets_sp1 = parse_miranda_output(seq_input_sp1, predicted_targets_sp1, 'conserved', shared_seeds, mature_sp1)
                targets_sp2 = parse_miranda_output(seq_input_sp2, predicted_targets_sp2, 'conserved', shared_seeds, mature_sp2)
            elif conservation == 'all_miRs':
                targets_sp1 = parse_miranda_output(seq_input_sp1, predicted_targets_sp1, 'all')
                targets_sp2 = parse_miranda_output(seq_input_sp2, predicted_targets_sp2, 'all')
        print('targets', len(targets_sp1), len(targets_sp2))
        
        # sort genes by UTR length
        UTR_length_sp1 = sort_genes_3UTR_length(UTR_sp1)
        UTR_length_sp2 = sort_genes_3UTR_length(UTR_sp2)
        print('UTR length', len(UTR_length_sp1), len(UTR_length_sp2))
        
        # create list of number of targets, normalized targets, sequence length for genes in CNV
        Sp1_CNV_sites, Sp2_CNV_sites, Sp1_CNV_norm_sites, Sp2_CNV_norm_sites, Sp1_CNV_length, Sp2_CNV_length = [], [], [], [], [], []
        
        # create list of number of targets, normalized targets, sequence length for genes in non CNV
        Sp1_nonCNV_sites, Sp2_nonCNV_sites, Sp1_nonCNV_norm_sites, Sp2_nonCNV_norm_sites, Sp1_nonCNV_length, Sp2_nonCNV_length = [], [], [], [], [], []        
        
        # loop over genes in targets sp1
        for gene in targets_sp1:
            # record only genes with short 3'UTR
            if UTR_length_sp1[gene] == 'short':
                # check if gene in CNV or not
                if CNV_status1[gene] == 'CNV':
                    Sp1_CNV_sites.append(targets_sp1[gene][0])
                    Sp1_CNV_norm_sites.append(targets_sp1[gene][2])
                    Sp1_CNV_length.append(targets_sp1[gene][1])
                elif CNV_status1[gene] == 'not_CNV':
                    Sp1_nonCNV_sites.append(targets_sp1[gene][0])
                    Sp1_nonCNV_norm_sites.append(targets_sp1[gene][2])
                    Sp1_nonCNV_length.append(targets_sp1[gene][1])
        
        # loop over genes in targets sp2
        for gene in targets_sp2:
            # record only genes with short 3'UTR
            if UTR_length_sp2[gene] == 'short':
                # check if gene in CNV or not
                if CNV_status2[gene] == 'CNV':
                    Sp2_CNV_sites.append(targets_sp2[gene][0])
                    Sp2_CNV_norm_sites.append(targets_sp2[gene][2])
                    Sp2_CNV_length.append(targets_sp2[gene][1])
                elif CNV_status2[gene] == 'not_CNV':
                    Sp2_nonCNV_sites.append(targets_sp2[gene][0])
                    Sp2_nonCNV_norm_sites.append(targets_sp2[gene][2])
                    Sp2_nonCNV_length.append(targets_sp2[gene][1])
                
                    
        print('CNV targets', len(Sp1_CNV_sites), len(Sp2_CNV_sites))
        print('non CNV targets', len(Sp1_nonCNV_sites), len(Sp2_nonCNV_sites))
        
        
        # write content to file
        newfile.write(species[i] + '\t' + species[j] + '\t')
        
        # write number of miRNAs if using shared miRNAs
        if conservation == 'shared_miRs':
            newfile.write(str(len(shared_seeds)) + '\t')
        
        newfile.write(str(len(Sp1_CNV_sites)) + '\t' + str(len(Sp2_CNV_sites)) + '\t')
        newfile.write(str(np.mean(Sp1_CNV_sites)) + '\t' + str(np.std(Sp1_CNV_sites) / math.sqrt(len(Sp1_CNV_sites))) + '\t')
        newfile.write(str(np.mean(Sp2_CNV_sites)) + '\t' + str(np.std(Sp2_CNV_sites) / math.sqrt(len(Sp2_CNV_sites))) + '\t')
        newfile.write(str(stats.ranksums(Sp1_CNV_sites, Sp2_CNV_sites)[1]) + '\t')
        newfile.write(str(len(Sp1_nonCNV_sites)) + '\t' + str(len(Sp2_nonCNV_sites)) + '\t')
        newfile.write(str(np.mean(Sp1_nonCNV_sites)) + '\t' + str(np.std(Sp1_nonCNV_sites) / math.sqrt(len(Sp1_nonCNV_sites))) + '\t')
        newfile.write(str(np.mean(Sp2_nonCNV_sites)) + '\t' + str(np.std(Sp2_nonCNV_sites) / math.sqrt(len(Sp2_nonCNV_sites))) + '\t')
        newfile.write(str(stats.ranksums(Sp1_nonCNV_sites, Sp2_nonCNV_sites)[1]) + '\t')
        newfile.write(str(np.mean(Sp1_CNV_norm_sites)) + '\t' + str(np.std(Sp1_CNV_norm_sites) / math.sqrt(len(Sp1_CNV_norm_sites))) + '\t')
        newfile.write(str(np.mean(Sp2_CNV_norm_sites)) + '\t' + str(np.std(Sp2_CNV_norm_sites) / math.sqrt(len(Sp2_CNV_norm_sites))) + '\t')
        newfile.write(str(stats.ranksums(Sp1_CNV_norm_sites, Sp2_CNV_norm_sites)[1]) + '\t')
        newfile.write(str(np.mean(Sp1_nonCNV_norm_sites)) + '\t' + str(np.std(Sp1_nonCNV_norm_sites) / math.sqrt(len(Sp1_nonCNV_norm_sites))) + '\t')
        newfile.write(str(np.mean(Sp2_nonCNV_norm_sites)) + '\t' + str(np.std(Sp2_nonCNV_norm_sites) / math.sqrt(len(Sp2_nonCNV_norm_sites))) + '\t')
        newfile.write(str(stats.ranksums(Sp1_nonCNV_norm_sites, Sp2_nonCNV_norm_sites)[1]) + '\t')              
        newfile.write(str(np.mean(Sp1_CNV_length)) + '\t' + str(np.std(Sp1_CNV_length) / math.sqrt(len(Sp1_CNV_length))) + '\t')
        newfile.write(str(np.mean(Sp2_CNV_length)) + '\t' + str(np.std(Sp2_CNV_length) / math.sqrt(len(Sp2_CNV_length))) + '\t')
        newfile.write(str(stats.ranksums(Sp1_CNV_length, Sp2_CNV_length)[1]) + '\t')
        newfile.write(str(np.mean(Sp1_nonCNV_length)) + '\t' + str(np.std(Sp1_nonCNV_length) / math.sqrt(len(Sp1_nonCNV_length))) + '\t')
        newfile.write(str(np.mean(Sp2_nonCNV_length)) + '\t' + str(np.std(Sp2_nonCNV_length) / math.sqrt(len(Sp2_nonCNV_length))) + '\t')
        newfile.write(str(stats.ranksums(Sp1_nonCNV_length, Sp2_nonCNV_length)[1]) + '\n')
        
        print('done writing regulation for {0} and {1}'.format(species[i], species[j]))

# close file after writing
newfile.close()



#########################



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
    
    
    
    
#############################



# -*- coding: utf-8 -*-
"""
Created on Wed Dec 23 14:44:01 2015
@author: RJovelin
"""

# create files to compare the number of mirna binding sites among species


from CNV_miRNAs import *
import os
import sys
import matplotlib.pyplot as plt

# usage python3 make_species_comp_files.py [options]
# [5UTR/CDS] use 5'UTR or CDS sites
# [True/False] use valid chromosomes (True) or all chromsomes (False)
# [long_CNVs/all_CNVs] use all CNVs or CNVs > 1Kb
# [CNV/not_CNV] compare CNV genes or non-CNV genes
# human_CNV_file choose the human CNV file

# get the region to consider to predict target sites [5UTr or CDS]
domain = sys.argv[1]
print(domain)

# get the option to keep genes on all chromos (False) or only on assembled 
# nuclear chromosomes only from the command
keep_valid_chromos = sys.argv[2]

# check if all chromos (including unplaced, unlocated, and MT) are used
# or if only valid chromos are used 
# note: only CNVs on valid chromos are reported in DGV, so if all chromos are
# used it may introduce a bias by calling non CNV genes genes that cannot be detected

if keep_valid_chromos == 'True':
    keep_valid_chromos = True
    chromos = 'valid_chromos'
elif keep_valid_chromos == 'False':
    keep_valid_chromos = False
    chromos = 'all_chromos'
print(keep_valid_chromos, chromos)


# get the type of CNVs to consider from the command
long_CNV = sys.argv[3]

if long_CNV == 'all_CNVs':
    cnv_length = 'CNV_all_length'
elif long_CNV == 'long_CNVs':
    cnv_length = 'CNV_greater_1Kb'
print(long_CNV, cnv_length)

# compare CNV genes or non-CNV genes
gene_type = sys.argv[4]

# get the human CNV file from the command, so that comparisons can be made 
# between all species and different version of the human DVG
human_CNV_file = sys.argv[5]

# get realease version of the DGV
release_version = human_CNV_file[human_CNV_file.index('GRCh37') : human_CNV_file.index('_CNV')]

# make a list of species names
species = ['H_sapiens', 'P_troglodytes', 'M_mulatta', 'M_musculus', 'B_taurus', 'G_gallus']

species_codes = {'H_sapiens': 'Hsa', 'P_troglodytes': 'Ptr', 'M_mulatta': 'Mml',
                 'M_musculus': 'Mmu', 'B_taurus': 'Bta', 'G_gallus':'Gga'}

# create a dict for human-species comparison with a list of list of sites for each species
# {Hsa_species_name : [[sites_Hsa], [sites_sp2]]}
species_data = {}

# loop over non-human species
for i in range(1, len(species)):
    # compare human to each other species
    print(species[0], species[i])
    # get genome files
    genome_sp1 = species[0] + '_genome.txt'
    genome_sp2 = species[i] + '_genome.txt'
    # get chromos files
    valid_chromos_sp1 = species[0] + '_valid_chromos.txt'
    valid_chromos_sp2 = species[i] + '_valid_chromos.txt'
    # get GFF annotation files
    GFF_sp1 = species[0] + '.gff3'
    GFF_sp2 = species[i] + '.gff3'
    print(genome_sp1, genome_sp2)
    print(valid_chromos_sp1, valid_chromos_sp2)
    print(GFF_sp1, GFF_sp2)

    # get targetscan sequence input file
    seq_input_sp1 = species[0] + '_' + domain + '_' + chromos + '_targetscan.txt'
    seq_input_sp2 = species[i] + '_' + domain + '_' + chromos + '_targetscan.txt'
    print(seq_input_sp1, seq_input_sp2)
        
    # get predictor output
    predicted_targets_sp1 = species[0] + '_' + domain + '_' + chromos + '_predicted_sites_targetscan.txt'
    predicted_targets_sp2 = species[i] + '_' + domain + '_' + chromos + '_predicted_sites_targetscan.txt'
    print(predicted_targets_sp1, predicted_targets_sp2)        
                
    # get UTR files
    UTR_sp1 = species[0] + '_3UTR_length_' + chromos + '.txt'      
    UTR_sp2 = species[i] + '_3UTR_length_' + chromos + '.txt'
    print(UTR_sp1, UTR_sp2)
         
    # get CNV files
    CNV_file_sp1 = human_CNV_file
    CNV_file_sp2 = species[i] + '_' + cnv_length + '_' + chromos + '.txt'
            
    # get mature mirna files
    mature_sp1 = species[0] + '_mature.txt'
    mature_sp2 = species[i] + '_mature.txt'        
    print(mature_sp1, mature_sp2)
        
    # get a set of shared seeds between species pair
    shared_seeds = find_conserved_mirna_families(mature_sp1, mature_sp2)
        
    # get CNV gene status
    CNV_status1 = sort_genes_CNV_status(CNV_file_sp1)
    CNV_status2 = sort_genes_CNV_status(CNV_file_sp2)
    print('CNV status', len(CNV_status1), len(CNV_status2))
        
    # parse predictor outputfile {gene: [N_sites, seq_length, N_sites/seq_length]}
    # consider only shared mirna families between species pairs
    # filter targets for shared mirnas
    targets_sp1 = parse_targetscan_output(seq_input_sp1, predicted_targets_sp1, 'conserved', shared_seeds, mature_sp1)
    targets_sp2 = parse_targetscan_output(seq_input_sp2, predicted_targets_sp2, 'conserved', shared_seeds, mature_sp2)
    print('targets', len(targets_sp1), len(targets_sp2))
        
    # sort genes by UTR length
    UTR_length_sp1 = sort_genes_3UTR_length(UTR_sp1)
    UTR_length_sp2 = sort_genes_3UTR_length(UTR_sp2)
    print('UTR length', len(UTR_length_sp1), len(UTR_length_sp2))
        
    # create list of normalized targets
    Sp1_sites, Sp2_sites =  [], []
        
    # loop over genes in targets sp1
    for gene in targets_sp1:
        # record only genes with short 3'UTR
        if UTR_length_sp1[gene] == 'short':
            # check if CNV or non-CNV genes should be recorded
            if CNV_status1[gene] == gene_type:
                # record normalized sites
                Sp1_sites.append(targets_sp1[gene][2])
                
    # loop over genes in targets sp2
    for gene in targets_sp2:
        # record only genes with short 3'UTR
        if UTR_length_sp2[gene] == 'short':
            # check if CNV or non-CNV genes should be recorded
            if CNV_status2[gene] == gene_type:
                Sp2_sites.append(targets_sp2[gene][2])
    
    print('targets', len(Sp1_sites), len(Sp2_sites))
    
    # populate dict
    comp_names = species_codes[species[0]] + '_' + species_codes[species[i]]
    print(comp_names)
    # copy lists to avoid modifying inner lists
    species_data[comp_names] = [Sp1_sites[:], Sp2_sites[:]]
    

for comp_names in species_data:
    print(comp_names, len(species_data[comp_names][0]), len(species_data[comp_names][1]))

# create a file to store the data
outputfile = 'Data_for_species_comparison_fig_' + domain + '_' + chromos + '_' + long_CNV + '_' + release_version + '_' + gene_type + '.txt' 
print(outputfile)
newfile = open(outputfile, 'w')
# write header
newfile.write('\t'.join(['Species1_Species2', 'Species', 'Sites']) + '\n')

# loop over dictionnary
for comp_names in species_data:
    # get species1 and species2
    species1, species2 = comp_names.split('_')
    # write comp_names, species_1 and corresponding sites
    newfile.write(comp_names + '\t' + species1 + '\t')
    newfile.write('\t'.join(list(map(lambda x: str(x), species_data[comp_names][0]))) + '\n')
    # write comp_names, species_2 and corresponding sites    
    newfile.write(comp_names + '\t' + species2 + '\t')
    newfile.write('\t'.join(list(map(lambda x: str(x), species_data[comp_names][1]))) + '\n')
    
# close file after writing
newfile.close()        
    
   














