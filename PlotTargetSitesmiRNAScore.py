# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 15:55:09 2016

@author: RJovelin
"""



# use this script to compare target sites between human CNV and non-CNv genes
# assigning different scores to miRNAs according to their expression level

# usage PlotTargetSitesmiRNAScore.py [options]
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


# create a dictionary with {species: {mirna accession : mature names pairs}}
AccessionNames = MatchmiRNAAccessionNumbers('name', miRBaseFile = 'miRNA.dat')
print('matches accessions with names')

# Record human accessions 
HumanAccessions = [accession for accession in AccessionNames['Homo_sapiens']]
print('human accessions', len(HumanAccessions))

# create a dictionary {mirna accession: expression level}
miRNAExpression = miRBAsemiRNAExpression('mirna_read_count.txt')
print('obtained mirna expression', len(miRNAExpression))

# remove non-human accessions 
to_delete = [mirna for mirna in miRNAExpression if mirna not in HumanAccessions]
for mirna in to_delete:
    del miRNAExpression[mirna]
print('removed non-human miRNAs', len(miRNAExpression))

# map mature names with expression level {mirna: expression}
MatureExpression = {}
for accession in miRNAExpression:
    for mature in AccessionNames['Homo_sapiens'][accession]:
        MatureExpression[mature] = miRNAExpression[accession]
print('matched mature names and expression level', len(MatureExpression))

# get the expression of human mirnas 
ExpressionLevel = [MatureExpression[mature] for mature in MatureExpression]
print('mirnas with expression', len(ExpressionLevel))
ExpressionLevel.sort()
print(len(ExpressionLevel), min(ExpressionLevel), max(ExpressionLevel), np.mean(ExpressionLevel), np.median(ExpressionLevel))

# create a histogram with expression level 
HistoCounts, HistoLimits = np.histogram(ExpressionLevel, range(0, 1000, 10))

# Create a score based on miRNA expression to weight the importance of mirna targets
# score is simply the bin position of the mirna expression
Score = {}
for mirna in MatureExpression:
    for i in range(0, len(HistoLimits) - 1):
        if MatureExpression[mirna] >= HistoLimits[i] and MatureExpression[mirna] < HistoLimits[i+1]:
            Score[mirna] = i+1
    if MatureExpression[mirna] >= HistoLimits[i+1]:
        Score[mirna] = len(HistoLimits)
assert len(Score) == len(MatureExpression), 'scores are not recorded for some mirnas'




# use this function to create dict from miranda outputs
def WeightTargetsMirandaOutput(targetscan_seq_input_file, predicted_targets, Score):
    '''
    (file, file, dict) -> dict
    Take the targetscan input sequence file, the miranda output, the dictionary with 
    mirna: expression score pairs and return a dictionary with gene as key and
    a list with the number of predicted target sites, sequence length and number of target sites
    normalized by sequence length for all miRNAs or for conserved miRNAs only.
    All target counts are weighted by a score to take into account the expression of the cognate miRNA
    '''
    
    # get the length of the sequences used to predict target sites {gene : seq_length}
    genes_length = get_domain_length_from_targetscan_input(targetscan_seq_input_file)
        
    # create a dict to store the target sites {gene: number of weighted targets} 
    target_counts = {}
    # create a dict with {gene : sequence length} pairs
    target_length = {}    
    
    # open file for reading
    infile = open(predicted_targets, 'r')
    # go through file
    for line in infile:
        if line.startswith('>>'):
            line = line.rstrip().split('\t')
            # get mirna, get rid of '>>' sign
            mirna = line[0][2:]
            # get target gene
            gene = line[1]
            # get sequence length
            seq_length = int(line[8])
            # populate dict with gene : sequence length pairs
            if gene not in target_length:
                target_length[gene] = seq_length
            # get positions
            positions = line[9].split()
            # populate dict
            if gene not in target_counts:
                # initialize value 
                target_counts[gene] = 0
            # count the number of target sites weighted by the score of the mirna
            target_counts[gene] += (len(positions) * Score[mirna])
    # close file after reading
    infile.close()
    
    # create a dict to store the number of targets {gene: [N_sites, seq_length, N_sites/seq_length]}
    targets = {}
    # loop over genes with predicted targets
    for gene in target_counts:
        # add number of target sites
        targets[gene] = [target_counts[gene]]
        # get the length of the sequence used to predict target sites
        seq_length = target_length[gene]
        # add sequence length to list
        targets[gene].append(seq_length)
        # add number of sites normalized by sequence length
        targets[gene].append(target_counts[gene] / seq_length)
    
    # count 0 for genes that do not have any targets but that have a region > mininum length
    # loop over genes in targetscan seq input
    for gene in genes_length:
        # check that sites are not already recorded
        if gene not in targets:
            # targets were not predicted
            # check that sequence is greater than 6 bp
            if genes_length[gene] >= 7:
                # populate dict
                targets[gene] = [0, genes_length[gene], 0]

    return targets


# get the seq input file
seq_input_file = 'H_sapiens_' + domain + '_' + chromos + '_targetscan.txt'
# get the predicted targets output file
predicted_targets = 'H_sapiens_' + domain + '_' + chromos + '_predicted_sites_miranda.txt'
# use DGV 2015 release 
CNV_file = 'H_sapiens_GRCh37_2015_CNV_all_length_valid_chromos.txt'


# record the number of miranda target sites for each gene weighted by the mirna expression score
#  {gene: [N_targets, Sequence_length, N_targets_normalized, CNV_status}}
Targets = WeightTargetsMirandaOutput(seq_input_file, predicted_targets, Score)

# get CNV gene status
CNV_status = sort_genes_CNV_status(CNV_file)
print('recorded CNV gene status')

# add CNV status
for gene in Targets:
    Targets[gene].append(CNV_status[gene])
    assert len(TargetsLowExp[gene]) == 4, 'gene in Targets does not have all required values'
print('added gene CNV status to each gene')



##### continue here





# create a dict with expression group as key and a list with the number of targets
# for CNv and non-CNV genes {expression: [[CNV], [non-CNV]]}
TargetsData = {}
# initialize dict
TargetsData['low'], TargetsData['moderate'], TargetsData['medium'], TargetsData['high'] = [[], []], [[], []], [[], []], [[], []]
for gene in TargetsLowExp:
    if TargetsLowExp[gene][-1] == 'CNV':
        # add number of normalized targets to cnv list 
        TargetsData['low'][0].append(TargetsLowExp[gene][2])
    elif TargetsLowExp[gene][-1] == 'not_CNV':
        # add number of normalized targets to non-cnv list
        TargetsData['low'][1].append(TargetsLowExp[gene][2])
for gene in TargetsModerateExp:
    if TargetsModerateExp[gene][-1] == 'CNV':
        # add number of normalized targets to cnv list
        TargetsData['moderate'][0].append(TargetsModerateExp[gene][2])
    elif TargetsModerateExp[gene][-1] == 'not_CNV':
        # add number of normalized targets to non-cnv list
        TargetsData['moderate'][1].append(TargetsModerateExp[gene][2])
for gene in TargetsMediumExp:
    if TargetsMediumExp[gene][-1] == 'CNV':
        # add number of normalized targets to cnv list
        TargetsData['medium'][0].append(TargetsMediumExp[gene][2])
    elif TargetsMediumExp[gene][-1] == 'not_CNV':
        # add numner of normnalized targets to non-cnv list
        TargetsData['medium'][1].append(TargetsMediumExp[gene][2])
for gene in TargetsHighExp:
    if TargetsHighExp[gene][-1] == 'CNV':
        # add number of normalized targets to cnv list
        TargetsData['high'][0].append(TargetsHighExp[gene][2])
    elif TargetsHighExp[gene][-1] == 'not_CNV':
        # add number of normalized targets to non-cnv list
        TargetsData['high'][1].append(TargetsHighExp[gene][2])
print('generated lists of target sites for CNV and non-CNV genes')


# perform stattistical tests between CNV and non-CNV genes
# create dicts to store results {expression group: P-value normalized targets}
CompTargets = {}
for group in TargetsData:
    Pval = stats.ranksums(TargetsData[group][0], TargetsData[group][1])[1]
    CompTargets[group] = Pval    
    print('{0} -  cnv: N = {1}, mean = {2}, non-cnv: N = {3}, mean = {4}, P = {5}'.format(group,
          len(TargetsData[group][0]), np.mean(TargetsData[group][0]), len(TargetsData[group][1]),
          np.mean(TargetsData[group][1]), Pval))    
print('compared CNV and non-CNV genes')

# make a list with all data
AllData = []
Groups = ['low', 'moderate', 'medium', 'high']
# loop over expression groups,  
for group in Groups:
    # add list of targets for cnv genes
    AllData.append(TargetsData[group][0])
    # add list of targets for non-cnv genes
    AllData.append(TargetsData[group][1])
print('data consolidated in array')


# create figure
fig = plt.figure(1, figsize = (4, 2.5))

# create subplot in figure
# add a plot to figure (N row, N column, plot N)
ax = fig.add_subplot(1, 1, 1)
# create a list of positions for the box plot    
BoxPositions = [0, 0.4, 0.9, 1.3, 1.8, 2.2, 2.7, 3.1]
# use a boxplot
bp = ax.boxplot(AllData, showmeans = True, showfliers = False, widths = 0.3,
                positions = BoxPositions, patch_artist = True) 

# color CNV and non-CNV boxes differently
CNVColor = ['#a6bddb','#74a9cf','#2b8cbe','#045a8d']
NonCNVColor = ['#99d8c9','#66c2a4','#2ca25f','#006d2c']
i, j = 0, 0    
# change box, whisker color to black
for box in bp['boxes']:
    # change line color
    box.set(color = 'black')
    if i % 2 == 0:
        # CNV data
        box.set(facecolor = CNVColor[j])
    else:
        box.set(facecolor = NonCNVColor[j])
    i += 1
    j = int(i / 2)
# change whisker color to black
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
    mean.set(marker = 'o', markeredgecolor = 'black', markerfacecolor = 'black', markersize = 3)
   

   
# set font for all text in figure
FigFont = {'fontname':'Arial'}   
    
# write label for y and x axis
ax.set_ylabel('Normalized number of miRNA\nsites per gene', color = 'black',  size = 8, ha = 'center', **FigFont)
ax.set_xlabel('miRNA expression level', color = 'black',  size = 8, ha = 'center', **FigFont)

# create list of labels and tick positions for the X axis
xtickpos = [0.2, 1.1, 2, 2.9]
ExGroups = ['Low', 'Moderate', 'Medium', 'High']
# write label for x axis
plt.xticks(xtickpos, ExGroups, ha = 'center', fontsize = 8, **FigFont)

## add a range for the Y and X axis
plt.ylim([0, 0.14])
plt.xlim([-0.25, 3.35])
    
# do not show lines around figure  
ax.spines["top"].set_visible(False)    
ax.spines["bottom"].set_visible(True)    
ax.spines["right"].set_visible(False)    
ax.spines["left"].set_visible(False)  
# offset the spines
for spine in ax.spines.values():
    spine.set_position(('outward', 5))
    
# add a light grey horizontal grid to the plot, semi-transparent, 
ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5, linewidth = 0.5)  
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
    labelsize = 8,
    direction = 'out') # ticks are outside the frame when bottom = 'on'  
      
# Set the tick labels font name
for label in ax.get_yticklabels():
    label.set_fontname('Arial')
   
# create a margin around the x axis
plt.margins(0.05)
   
# annotate Graph with significance level
Pvalues = []
for group in Groups:
    # get the significance level for target sites
    if CompTargets[group] >= 0.05:
        Pvalues.append('')
    elif CompTargets[group] < 0.05 and CompTargets[group] >= 0.01:
        Pvalues.append('*')
    elif CompTargets[group] < 0.01 and CompTargets[group] >= 0.001:
        Pvalues.append('**')
    elif CompTargets[group] < 0.001:
        Pvalues.append('***')

# create list of Y and X positions to annotate figure with significance level
if domain == '3UTR':
    # make a list of Y positions
    Ypos = [0.13, 0.10, 0.09, 0.09]
    Xpos = [0.2, 1.1, 2, 2.9]


# annotate figure with significance levels
for i in range(len(Pvalues)):
    ax.text(Xpos[i], Ypos[i], Pvalues[i], horizontalalignment = 'center', verticalalignment = 'center', color = 'black', size = 8)


# add legend relative to ax1 using ax1 coordinates
C = mpatches.Patch(facecolor = '#a6bddb', edgecolor = 'black', linewidth = 1, label= 'CNV')
N = mpatches.Patch(facecolor = '#99d8c9', edgecolor = 'black', linewidth = 1, label= 'non-CNV')
ax.legend(handles = [C, N], loc = (0.2, 1), fontsize = 8, frameon = False, ncol = 2)

## build outputfile with arguments
#outputfile = 'truc_' + domain + '_' + chromos + '_' + cnv_length
#print(outputfile)

# save figure
fig.savefig('truc.pdf', bbox_inches = 'tight')

