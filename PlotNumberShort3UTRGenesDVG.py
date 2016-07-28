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

# sort genes based on 3' UTR length
UTR_length = sort_genes_3UTR_length(UTR_file, L)
print(len(UTR_length))

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

# create a dict {release: {study: {set of cnv genes}}}
StudiesCNVGenes = {}

# make a dictionary to match each study to a pubmed ID in each release
# {release: {reference: pubmedid}}
References = {}
for filename in CNV_files:
    # get release version
    release_version = filename[filename.index('GRCh37'): filename.index('_CNV')]
    References[release_version] = {}
    ref = get_DGV_references(filename)
    References[release_version] = dict(ref)
print('got references')

# get synonym names for all genes {gene name : [list of synonyms]}
synonyms = get_synonyms('H_sapiens.gff3')
print('got synonymous names', len(synonyms))

# get the CDS sequences of the longest mRNAs for each gene {gene : sequence}
CDS_seq = extract_CDS_sequences('H_sapiens.gff3', 'H_sapiens_genome.txt', 'H_sapiens_valid_chromos.txt', keep_valid_chromos)
print('extracted CDS sequences', len(CDS_seq))  


# get the CNV genes for each study of each release of the DGV
# loop over CNV files
for filename in CNV_files:
    print(filename)
    # sort genes based on CNV status
    CNV_status = sort_genes_CNV_status(filename)
    print(len(CNV_status))
    # get release version
    release_version = filename[filename.index('GRCh37'): filename.index('_CNV')]
    print(release_version)
    # initialize outer dict
    StudiesCNVGenes[release_version] = {}
    # get the set of CNV genes for each study of that release
    for study in References[release_version]:
        CNV_genes = get_human_CNV_genes_single_study(filename, study, 'all')
        StudiesCNVGenes[release_version][study] = set(CNV_genes)
print('got CNV genes for each study')        
        
    
# create a dict {release: {study: {gene: CNV status}}}    
CNV_status = {}    
for release in StudiesCNVGenes:
    # initialize outer dict
    CNV_status[release] = {}
    for study in StudiesCNVGenes[release]:
        # initialize dict
        StudiesCNVGenes[release][study] = {}
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
    











# make a dict with the number of CNV genes and non-CNV genes for each study
CNVNum = {}
for study in CNVTargets:
    CNVNum[study] = [0, 0]
    for gene in CNVTargets[study]:
        if CNVTargets[study][gene][-1] == 'CNV':
            CNVNum[study][0] += 1
        elif CNVTargets[study][gene][-1] == 'not_CNV':
            CNVNum[study][1] += 1
print('got CNV gene counts for each study')
for study in CNVNum:
    print(study, CNVNum[study][0], CNVNum[study][1])

# create a dict of {species name : [[CNV], [non-CNV]]}
CNVData = {} 
# loop over studies, get the number of normalized sites for CNV and non-CNV genes
for study in CNVTargets:
    # initialise list value
    CNVData[study] = [[], []]
    # populate inner lists with number of miRNA target sites per nucleotide
    for gene in CNVTargets[study]:
        if CNVTargets[study][gene][-1] == 'CNV':
            CNVData[study][0].append(CNVTargets[study][gene][2])
        elif CNVTargets[study][gene][-1] == 'not_CNV':
            CNVData[study][1].append(CNVTargets[study][gene][2])
print('generated lists of target sites for CNV and non-CNV genes')


# make a list of all data for each study
AllData = []
# make a list of species names to loop from
StudyNames = ['Suktitipat_et_al_2014', 'Alsmadi_et_al_2014', 'John_et_al_2014', 'Thareja_et_al_2015']
Populations = ['$Thai^a$', '$Kuwaiti^b$', '$Kuwaiti^c$', '$Kuwaiti^d$']

# loop over study in studies list and populate data lists, keeping the same order for targetscan and miranda
for study in StudyNames:
    # append list of target sites for CNV genes
    AllData.append(CNVData[study][0])
    # append list of target sites for non-CNV genes
    AllData.append(CNVData[study][1])
print('data consolidated in array')


# boostrap CNV and non-CNV genes to compare miRNA targets
# create a dictionary {study: {CNV_status: {num: gene}}} 
ToSampleFrom = {}
# loop over studies
for study in CNVTargets:
    # initialize dict
    ToSampleFrom[study] = {}
    ToSampleFrom[study]['CNV'], ToSampleFrom[study]['not_CNV'] = {}, {}
    # initialize counters for CNV and non-CNV genes
    i, j = 0, 0
    # loop over genes for that study
    for gene in CNVTargets[study]:
        if CNVTargets[study][gene][-1] == 'CNV':
            # populate dict with num : gene pair and update counter
            ToSampleFrom[study]['CNV'][i] = gene
            i += 1
        elif CNVTargets[study][gene][-1] == 'not_CNV':
            # populate dict with num : gene pair and update counter
            ToSampleFrom[study]['not_CNV'][j] = gene
            j += 1
              
# check that all genes have been assigned to a number
for study in ToSampleFrom:
    assert len(ToSampleFrom[study]['CNV']) == CNVNum[study][0], 'CNV genes counts do not match'
    assert len(ToSampleFrom[study]['not_CNV']) == CNVNum[study][1], 'non-CNV genes counts do not match'


# create a dict for each study with a list with numbers of each different outcomes when comparing targets in CNV and non-CNV genes
# {study: [# replicates CNV > non-CNV, # replicates CNV < non-CNV, # replicates no differences]}
BootStrap = {}
# initialize list values
for study in ToSampleFrom:
    BootStrap[study] = [0, 0, 0]

# loop over studies in dict to sample from
for study in ToSampleFrom:
    print('bootstraping', study)
    # set number of replicates
    replicates = 10000
    while replicates != 0:
        # make list of targets for CNV and non-CNV genes
        repCNVtargets, repNonCNVtargets = [], []
        # draw 800 CNV genes and 800 non-CNV genes with replacement
        for i in range(800):
            # draw a random CNV gene
            j = random.randint(0, len(ToSampleFrom[study]['CNV']) - 1)
            k = random.randint(0, len(ToSampleFrom[study]['not_CNV']) - 1)
            # get the corresponding genes
            gene1 = ToSampleFrom[study]['CNV'][j]
            gene2 = ToSampleFrom[study]['not_CNV'][k]            
            # get the the number of targets for these 2 genes
            assert CNVTargets[study][gene1][-1] == 'CNV', 'random gene should be CNV'
            assert CNVTargets[study][gene2][-1] == 'not_CNV', 'random gene should be non-CNV'
            repCNVtargets.append(CNVTargets[study][gene1][2])
            repNonCNVtargets.append(CNVTargets[study][gene2][2])
        # make sure that the correct numbers of genes is drawn
        assert len(repCNVtargets) == 800, '800 CNV genes should be drawn'
        assert len(repNonCNVtargets) == 800, '800 non-CNV genes should be drawn'
        # compare CNV and non-CNV genes
        Pval = stats.ranksums(repCNVtargets, repNonCNVtargets)[1]
        # check significance
        if Pval >= 0.05:
            # difference is not significance
            BootStrap[study][2] += 1
        elif Pval < 0.05:
            # difference is significance, check if CNV genes have a greater number of targets
            if np.mean(repCNVtargets) > np.mean(repNonCNVtargets):
                BootStrap[study][0] += 1
            elif np.mean(repCNVtargets) < np.mean(repNonCNVtargets):
                BootStrap[study][1] += 1
        # update replicate number
        replicates -= 1
print('done with boostraping')        

      
# create parallel list of proportions 
Greater, Lower, Nodiff = [], [] ,[]
for study in StudyNames:
    Greater.append(BootStrap[study][0] / sum(BootStrap[study]))
    Lower.append(BootStrap[study][1] / sum(BootStrap[study]))
    Nodiff.append(BootStrap[study][2] / sum(BootStrap[study]))

# create a list with all the proportion lists
Proportions = [Greater, Lower, Nodiff]

# print results to screen
for study in BootStrap:
    print(BootStrap[study])



# create figure
fig = plt.figure(1, figsize = (8, 3))

# create a function to format the subplots
def CreateAx(Columns, Rows, Position, Data, figure, Title, YMax, LabelNames, XScale, GraphType):
    '''
    (int, int, int, list, figure_object, str, int, list, list)
    Take the number of a column, and rows in the figure object and the position of
    the ax in figure, a list of data, a title, a maximum value for the Y axis,
    a list with species names and list of X axis tick positions and return an
    ax instance in the figure
    '''    
    
    # create subplot in figure
    # add a plot to figure (N row, N column, plot N)
    ax = figure.add_subplot(Rows, Columns, Position)
    
    if GraphType == 'box':
        # create a list of positions for the box plot    
        BoxPositions = [0, 0.4, 0.9, 1.3, 1.8, 2.2, 2.7, 3.1]
        # use a boxplot
        bp = ax.boxplot(Data, showmeans = True, showfliers = False, widths = 0.3,
                        positions = BoxPositions, patch_artist = True) 
    
        # color CNV and non-CNV boxes differently
        i = 0    
        # change box, whisker color to black
        for box in bp['boxes']:
            # change line color
            box.set(color = 'black')
            if i % 2 == 0:
                # CNV data, color box in grey
                box.set(facecolor = '#a6cee3')
            else:
                box.set(facecolor = '#b2df8a')
            i += 1
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
    
        # add legend
        C = mpatches.Patch(facecolor = '#a6cee3' , edgecolor = 'black', linewidth = 1, label= 'CNV')
        N = mpatches.Patch(facecolor = '#b2df8a' , edgecolor = 'black', linewidth = 1, label= 'non-CNV')
        plt.legend(handles = [C, N], loc = (0, 1), fontsize = 8, frameon = False, ncol = 2)
    
  
    elif GraphType == 'bar':
        # get the list of proportions 
        greater, lower, nodiff = Data[0], Data[1], Data[2]
        # make a list with added values for nodiff and greater
        added = []
        for i in range(len(greater)):
            added.append(nodiff[i] + greater[i])
        # Create a bar plot for proportions of replicates with CNV no diff on top of CNV lower
        ax.bar([0, 0.4, 0.8, 1.2], nodiff, width = 0.3, label = 'No difference', color= '#f7f7f7')
        # Create a bar plot for proportions of replicates with CNV greater on top of no diff
        ax.bar([0, 0.4, 0.8, 1.2], greater, width = 0.3, bottom = nodiff, label = 'CNV > non-CNV', color= '#ef8a62')
        # Create a bar plot for proportions of replicates with CNV lower on top of CNV greater
        ax.bar([0, 0.4, 0.8, 1.2], lower, width = 0.3, bottom= added, label = 'CNV < non-CNV', color = '#67a9cf')
 
        # add legend
        N = mpatches.Patch(facecolor = '#f7f7f7' , edgecolor = 'black', linewidth = 1, label= 'No diff.')
        G = mpatches.Patch(facecolor = '#ef8a62' , edgecolor = 'black', linewidth = 1, label= 'CNV greater')
        L = mpatches.Patch(facecolor = '#67a9cf' , edgecolor = 'black', linewidth = 1, label= 'CNV lower')
        plt.legend(handles = [N, G, L], loc = (0, 1), fontsize = 8, frameon = False, ncol = 3)

    # write title   
    ax.set_title(Title + '\n\n', size = 8)
    
    # set font for all text in figure
    FigFont = {'fontname':'Arial'}   
    
    # write label for y axis
    if GraphType == 'box':
        ax.set_ylabel('Normalized number of miRNA\nsites per gene', color = 'black',  size = 8, ha = 'center', **FigFont)
    elif GraphType == 'bar':
        ax.set_ylabel('Proportion of replicates', color = 'black', size = 8, ha = 'center', **FigFont)

    # write label for x axis
    plt.xticks(XScale, LabelNames, ha = 'center', fontsize = 8, **FigFont)

    # add a range for the Y amd X axes
    if GraphType == 'box':
        plt.ylim([0, YMax])
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
    
    return ax      



# plot boxplots for predictor in 1st subplot
# get title based on predictor algorithm
if predictor == 'targetscan':
    figtitle = 'TargetScan'
elif predictor == 'miranda':
    figtitle = 'miRanda'
ax1 = CreateAx(2, 1, 1, AllData, fig, figtitle, 0.45, Populations, [0.2, 1.1, 2, 2.9], 'box')
# plot bars in 2nd subplot
ax2 = CreateAx(2, 1, 2, Proportions, fig, figtitle, 0.45, Populations, [0.15, 0.55, 0.95, 1.35], 'bar')


# annotate Graph with significance level
Pvalues = []
for study in StudyNames:
    # compare mirna targets between CNV and non-CNV genes
    P = stats.ranksums(CNVData[study][0], CNVData[study][1])[1]
    if P >= 0.05:
        Pvalues.append('')
    elif P < 0.05 and P >= 0.01:
        Pvalues.append('*')
    elif P < 0.01 and P >= 0.001:
        Pvalues.append('**')
    elif P < 0.001:
        Pvalues.append('***')
# create list of Y and X positions to annotate figure with significance level
Ypos = [0.42, 0.42, 0.42, 0.42]
Xpos = [0.2, 1.1, 2, 2.9]
for i in range(len(Pvalues)):
    ax1.text(Xpos[i], Ypos[i], Pvalues[i], horizontalalignment = 'center',
             verticalalignment = 'center', color = 'black', size = 8)


# add subplot label
ax1.text(-1, 0.48, 'A', horizontalalignment = 'center',
         verticalalignment = 'center', color = 'black', size = 10)
ax1.text(3.5, 0.48, 'B', horizontalalignment = 'center',
         verticalalignment = 'center', color = 'black', size = 10)


# make sure subplots do not overlap
plt.tight_layout()

# get outputfile
outputfile = 'PlotSinglePops_' + predictor + '_' + domain + '_' + chromos + '_' + cnv_length 
print(outputfile)

# save figure
fig.savefig(outputfile + '.eps', bbox_inches = 'tight')
       
Contact GitHub API Training Shop Blog About
© 2016 GitHub, Inc. Terms Privacy Security Status Help








$$$$$$$$$$$$$$$





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
