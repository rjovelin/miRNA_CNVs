# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 16:53:07 2016

@author: RJovelin
"""


# use this script to plot the frequency of studies for which the number of
# miRNA sites for CNV genes is greater, lower or similar to non-CNV genes

# usage PlotFreqStudies.py [options]
# [3UTR/5UTR/CDS]: gene domain to consider


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


# get the region to consider to predict target sites [3UTR or 5UTr or CDS]
domain = sys.argv[1]
print(domain)

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


# make a dictionary of species names : species code
species_codes = {'H_sapiens': 'Hsa', 'P_troglodytes': 'Ptr', 'M_mulatta': 'Mml',
                 'M_musculus': 'Mmu', 'B_taurus': 'Bta', 'G_gallus':'Gga'}



















######################


# usage single_study_average_target_sites.py [parameters]
# [3UTR/5UTR/CDS] choose the region to analyse
# [targetscan/miranda] choose the target predictor
# [True/False] use valid chromos (True) or all chromos
# [long_CNVs/all_CNVs] use all CNVs or CNVs > 1 Kb
# CNV_file Choose DGV CNV file
# minimum_cnv minimum number of cnv genes in single study


from CNV_miRNAs import *
import os
import sys
import numpy as np
from scipy import stats
import math

# get the region to consider to predict target sites [3UTR or 5UTr or CDS]
domain = sys.argv[1]
print(domain)

# get predictor from command [targetscan or miranda]
predictor = sys.argv[2]
print(predictor)

# get the option to keep genes on all chromos (False) or only on assembled 
# nuclear chromosomes only (True) from the command
keep_valid_chromos = sys.argv[3]
if keep_valid_chromos == 'True':
    keep_valid_chromos = True
    chromos = 'valid_chromos'
elif keep_valid_chromos == 'False':
    keep_valid_chromos = False
    chromos = 'all_chromos'
print(keep_valid_chromos, chromos)


# get the option to call a CNV if CNV length > 1 Kb (long_CNVs)
# or to include all CNVs regardless of length (all_CNVs)
long_CNV = sys.argv[4]
if long_CNV == 'all_CNVs':
    cnv_length = 'CNV_all_length'
    CNV_size = 'all'
elif long_CNV == 'long_CNVs':
    cnv_length = 'CNV_greater_1Kb'
    CNV_size = 'long'
print(long_CNV, cnv_length, CNV_size)

# get the CNV file from the command line
# ie. all_CNV or CNVs > 1 Kb can be used for any DGV release
CNV_file = sys.argv[5]
print(CNV_file)

# get the minimum number of cnv genes that a single study must have
minimum_cnv = int(sys.argv[6])
print(minimum_cnv)


# check if all chromos (including unplaced, unlocated, and MT) are used
# or if only valid chromos are used 
# note: only CNVs on valid chromos are reported in DGV, so if all chromos are
# used it may introduce a bias by calling non CNV genes genes that cannot be detected

# get the file with 3'UTR length
UTR_file = 'H_sapiens_3UTR_length_' + chromos + '.txt'
print(UTR_file)
# get targetscan sequence input file
targetscan_seq_input_file = 'H_sapiens_' + domain + '_' + chromos + '_targetscan.txt'
print(targetscan_seq_input_file)

# get the outputfile with predicted target sites
predicted_target_file = 'H_sapiens_' + domain + '_' + chromos + '_predicted_sites_' + predictor + '.txt'    
print(predicted_target_file)
# make a dictionary with {gene :[targets, seq_length, normalized_targets]}

# check predictor
if predictor == 'targetscan':
    predicted_targets = parse_targetscan_output(targetscan_seq_input_file, predicted_target_file, 'all')
elif predictor == 'miranda':
    predicted_targets = parse_miranda_output(targetscan_seq_input_file, predicted_target_file, 'all')

# get release version
if '2013' in CNV_file:
    release_version = 'GRCh37_' + CNV_file[CNV_file.rindex('_') + 1: -7]
else:
    release_version = 'GRCh37_' + CNV_file[CNV_file.rindex('_') + 1: -10]
print(release_version)

# get outputfile
outputfile = 'H_sapiens_single_study_targets_' + domain + '_' + cnv_length + '_' + chromos + '_' + predictor + '_' + release_version + '.txt'
print(outputfile)

# open file for writing
newfile = open(outputfile, 'w')

# write header to file
newfile.write('\t'.join(['Study', 'N_CNV_genes', 'CNV_mean_targets', 'CNV_SEM_targets',  
                        'N_nonCNV_genes', 'nonCNV_mean_targets', 'nonCNV_SEM_targets',
                        'P_diff_targets', 'CNV_mean_normalized_targets', 'CNV_SEM_normalized_targets',
                        'nonCNV_mean_normalized_targets', 'nonCNV_SEM_normalized_targets', 'P_diff_normalized_targets',
                        'CNV_mean_seq_length', 'CNV_SEM_seq_length', 'nonCNV_mean_seq_length',
                        'nonCNV_SEM_seq_length', 'P_seq_length', 'Spearman_rho_targets_X_length', 'P_targets_X_length']) + '\n')                        

# get the dictionaries of {reference: pubmedid} for studies reported in the CGV
references = get_DGV_references(CNV_file)

# create a lamda function to transform value into string
Gstr = lambda x: str(x)

# loop over study
for study in references:
    print(study)
    # get the set of CNV genes corresponding to that study
    CNV_genes = get_human_CNV_genes_single_study(CNV_file, study, CNV_size)
    print('# CNV genes', len(CNV_genes)) 
    
    # check that study includes minimum number of cnv genes
    if len(CNV_genes) >= minimum_cnv:
        # make temporary cnv_file with CNV genes extracted from DGV
        tempfile = open('Temp_cnv_file.txt', 'w')
        # dump all CNV genes
        for gene in CNV_genes:
            tempfile.write(gene + '\n')
        tempfile.close()
    
        # get CNV gene status
        CNV_status = get_genes_CNV_status('H_sapiens.gff3', 'H_sapiens_genome.txt', 'H_sapiens_valid_chromos.txt', keep_valid_chromos, 'Temp_cnv_file.txt')    
        print('genes with CNV status', len(CNV_status))
    
        # make a temporary file with CNV status of all genes     
        tempfile = open('Temp_CNV_status_file.txt', 'w')
        tempfile.write('gene\tCNV_status\n')
        for gene in CNV_status:
            tempfile.write(gene + '\t' + CNV_status[gene] + '\n')
        tempfile.close()
    
        # make temp summary file with targets and cnv status
        make_summary_table_target_sites(predicted_targets, 'Temp_CNV_status_file.txt', 'Temp_summary_targets.txt')
                
        # count the number of CNV genes for that study
        Num_cnv_genes = 0
        # open temp summary file for reading
        infile = open('Temp_summary_targets.txt', 'r')
        # skip header
        infile.readline()
        # loop over file
        for line in infile:
            line = line.rstrip()
            if line != '':
                line = line.split()
                if line [-1] == 'CNV':
                    Num_cnv_genes += 1
        #close file after reading
        infile.close()
        print('Num CNV genes', Num_cnv_genes)
    
        # check that study includes minimum number of cnv genes
        if Num_cnv_genes >= minimum_cnv:
            # parse the summary table into a list
            regulation = compare_miRNA_regulation('Temp_summary_targets.txt')
    
            # write regulation to file
            newfile.write(study + '\t')
            newfile.write('\t'.join(list(map(Gstr, regulation))) + '\n')
    
            print('done writing regulation for {0}'.format(study))

# close file after writing
newfile.close()






















#################



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
    