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
# minimum_cnv =  minimum number of cnv genes in single study
MinimumCNVGenes = 500

Releases = ['GRCh37_2013-05', 'GRCh37_2013-07', 'GRCh37_2014', 'GRCh37_2015']
labelnames = ['2013a', '2013b', '2014', '2015']



# get the number of target sites for each study of each version of the DGV

# get UTR file
UTR_file = 'H_sapiens_3UTR_length_' + chromos + '.txt'
print(UTR_file)

# get synonym names for all genes {gene name : [list of synonyms]}
synonyms = get_synonyms('H_sapiens.gff3')
print('got synonymous names', len(synonyms))

# get the CDS sequences of the longest mRNAs for each gene {gene : sequence}
CDS_seq = extract_CDS_sequences('H_sapiens.gff3', 'H_sapiens_genome.txt', 'H_sapiens_valid_chromos.txt', keep_valid_chromos)
print('extracted CDS sequences', len(CDS_seq))  

# make a list of DGV files
DGVFiles = ['GRCh37_hg19_variants_2013-05-31.txt', 'GRCh37_hg19_variants_2013-07-23.txt',
            'GRCh37_hg19_variants_2014-10-16.txt', 'GRCh37_hg19_variants_2015-07-23.txt']

# get the CNV genes for each study of each release of the DGV
# create a dict {release: {study: {set of cnv genes}}}
StudiesCNVGenes = {}
# loop over DGV files
for filename in DGVFiles:
    print(filename)
    # get release version
    if '2013' in filename:
        release_version = filename[:filename.index('_hg19')] + '_' + filename[filename.index('variants_') + len('variants_'): -7]
    else:
        release_version = filename[:filename.index('_hg19')] + '_' + filename[filename.index('variants_') + len('variants_'): -10]
    print(release_version)
    # initialize outer dict
    StudiesCNVGenes[release_version] = {}
    # get the set of CNV genes for each study of that release
    for study in References[release_version]:
        CNV_genes = get_human_CNV_genes_single_study(filename, study, 'all')
        StudiesCNVGenes[release_version][study] = set(CNV_genes)
print('got CNV genes for each study')        
        
# remove studies if the number of CNV genes < MinimumCNVGenes
for release in StudiesCNVGenes:
    to_remove = []
    for study in StudiesCNVGenes[release]:
        if len(StudiesCNVGenes[release][study]) < MinimumCNVGenes:
            to_remove.append(study)
    if len(to_remove) != 0:
        for study in to_remove:
            del StudiesCNVGenes[release][study]
# remove version if all studies were removed
to_remove = []
for release in StudiesCNVGenes:
    if len(StudiesCNVGenes[release]) == 0:
        to_remove.append(release)
if len(to_remove) != 0:
    for release in to_remove:
        del StudiesCNVGenes[release]
# check that all 4 releases are kept
assert len(StudiesCNVGenes) == 4, 'not all releases are recorded'


# get the CNV status of all genes used to predict target sites for each study of each release
# create a dict {release: {study: {gene: CNV status}}}    
CNV_status = {}    
for release in StudiesCNVGenes:
    # initialize outer dict
    CNV_status[release] = {}
    for study in StudiesCNVGenes[release]:
        # initialize dict
        CNV_status[release][study] = {}
        # loop over genes in CDS seq
        for gene in CDS_seq:
            # set boolean
            is_cnv = False
            # ask if gene in CNV gene
            if gene in StudiesCNVGenes[release][study] or gene.upper() in StudiesCNVGenes[release][study]:
                CNV_status[release][study][gene] = 'CNV'
            else:
                # ask if any of the gene synonyms are in CNVs
                for name in synonyms[gene]:
                    # check if name in CNV
                    if name in StudiesCNVGenes[release][study] or name.upper() in StudiesCNVGenes[release][study]:
                        # update boolean
                        is_cnv = True
                # check if gene in CNV
                if is_cnv == True:
                    CNV_status[release][study][gene] = 'CNV'
                elif is_cnv == False:
                    CNV_status[release][study][gene] = 'not_CNV'
print('sorted genes according to CNV status')


for release in CNV_status:
    for study in CNV_status[release]:
        print(release, study, len(CNV_status[release][study]))




# get the number of target sites for each gene
# get targetscan sequence input file
targetscan_seq_input_file = 'H_sapiens_' + domain + '_' + chromos + '_targetscan.txt'
print(targetscan_seq_input_file)
# get the outputfile with predicted target sites
TargetScanTargetsFile = 'H_sapiens_' + domain + '_' + chromos + '_predicted_sites_targetscan.txt'   
MirandaTargetsFile = 'H_sapiens_' + domain + '_' + chromos + '_predicted_sites_miranda.txt'   
# make a dictionary with {gene :[targets, seq_length, normalized_targets]}
TargetsTargetscan = parse_targetscan_output(targetscan_seq_input_file, TargetScanTargetsFile, 'all')
TargetsMiranda = parse_miranda_output(targetscan_seq_input_file, MirandaTargetsFile, 'all')


# create dicts for targetscan and miranda predictors with target sites and CNV status
# {release: {study: {gene: [targets, seq_length, normalized_targets, CNV_status]}}}
CNVTargetsTargetscan, CNVTargetsMiranda = {}, {}
# loop over release in CNV status
for release in CNV_status:
    # initialize inner dict
    CNVTargetsTargetscan[release], CNVTargetsMiranda[release] = {}, {}
    for study in CNV_status[release]:
        # initialize inner dict
        CNVTargetsTargetscan[release][study], CNVTargetsMiranda[release][study] = {}, {}
        # loop over genes, populate dict combining targets and CNv status 
        for gene in CNV_status[release][study]:
            # get CNV status and add it to the list of targets, without modifying the original list
            status = CNV_status[release][study][gene]
            CNVTargetsTargetscan[release][study][gene] = list(TargetsTargetscan[gene]) + [status]
            CNVTargetsMiranda[release][study][gene] = list(TargetsMiranda[gene]) + [status]
print('combined CNV status and miRNA targets for each gene in each study')            


# compare the normalized number of target sites between CNV and non-CNV genes for each study and DGV version
# count the number of studies with CNV genes having greater, lowe or similar number targets as non-CNV genes
# {release: [CNV_greater, CNV_lower, NO_diff]}
CompTargetscan, CompMiranda = {}, {}
for release in CNVTargetsTargetscan:
    # initialize dict
    CompTargetscan[release] = [0, 0, 0]
    # loop over studies
    for study in CNVTargetsTargetscan[release]:
        # make lists of target sites for CNV amd non-CNV genes
        cnvtargets = [CNVTargetsTargetscan[release][study][gene][2] for gene in CNVTargetsTargetscan[release][study] if CNVTargetsTargetscan[release][study][gene][-1] == 'CNV']
        noncnvtargets = [CNVTargetsTargetscan[release][study][gene][2] for gene in CNVTargetsTargetscan[release][study] if CNVTargetsTargetscan[release][study][gene][-1] == 'not_CNV']
        # compare the number of target sites
        P = stats.ranksums(cnvtargets, noncnvtargets)[1]
        # update counters
        if P < 0.05:
            # compare means
            if np.means(cnvtargets) > np.means(noncnvtargets):
                CompTargetscan[release][0] += 1
            elif np.means(cnvtargets) < np.means(noncnvtargets):
                CompTargetscan[release][1] += 1
        elif P >= 0.05:
            CompTargetscan[release][2] += 1
for release in CNVTargetsMiranda:
    # initialize dict
    CompMiranda[release] = [0, 0, 0]
    # loop over studies
    for study in CNVTargetsMiranda[release]:
        # make lists of target sites for CNV and non-CNV genes
        cnvtargets = [CNVTargetsMiranda[release][study][gene][2] for gene in CNVTargetsMiranda[release][study] if CNVTargetsMiranda[release][study][gene][-1] == 'CNV']
        noncnvtargets = [CNVTargetsMiranda[release][study][gene][2] for gene in CNVTargetsMiranda[release][study] if CNVTargetsMiranda[release][study][gene][-1] == 'non_CNV']
        # compare the number of target sites
        P = stats.ranksums(cnvtargets, noncnvtargets)[1]
        # update counters
        if P < 0.05:
            # compare means
            if np.means(cnvtargets) > np.means(noncnvtargets):
                CompMiranda[release][0] += 1
            elif np.means(cnvtargets) < np.means(noncnvtargets):
                CompMiranda[release][1] += 1
        elif P >= 0.05:
            CompMiranda[release][2] += 1
print('compared mean target sites between CNV and non-CNV genes')        






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
    