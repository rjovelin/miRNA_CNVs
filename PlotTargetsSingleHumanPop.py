# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 14:10:59 2016

@author: RJovelin
"""


# use this script to compare the number target sites per nucleotide between CNV and non-CNV genes
# for single human populations 

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


# usage PlotTargetsSingleHumanPop.py [options]
# [3UTR/5UTR/CDS] choose the region to analyse
# -[targetscan/miranda] : the predictor algorithm used to predict target sites


# get the region to consider to predict target sites [3UTR or 5UTr or CDS]
domain = sys.argv[1]
print(domain)
# get predictor from command [targetscan or miranda]
predictor = sys.argv[2]
print(predictor)

# keep genes on assembled chromsomes
keep_valid_chromos = True
chromos = 'valid_chromos'
# use all CNVs
cnv_length = 'CNV_all_length'

# get the CNV files with single human population references
CNV_file = 'GRCh37_hg19_variants_2016-05-15.txt'

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
if predictor == 'targetscan':
    predicted_targets = parse_targetscan_output(targetscan_seq_input_file, predicted_target_file, 'all')
elif predictor == 'miranda':
    predicted_targets = parse_miranda_output(targetscan_seq_input_file, predicted_target_file, 'all')
print('# targets', len(predicted_targets))


# make a dictionary {reference: pubmedid} for studies reported in the CGV
references = get_DGV_references(CNV_file)

# make a list of pubmed IDs for single population studies
pubmed = [['Suktitipat_et_al_2014', '25118596'], ['John_et_al_2014', '26484159'], ['Thareja_et_al_2015', '25765185'], ['Alsmadi_et_al_2014', '24896259']]

# check that reference names correspond to the pubmed IDs
infile = open(CNV_file)
for line in infile:
    for i in range(len(pubmed)):
        if pubmed[i][0] in line:
            assert pubmed[i][1] in line, 'study ref does not correspond to expected pubmed ID'
infile.close()


# create a set of CNV genes for each study {study: {gene1, gene2...}}
StudiesCNV = {}
for i in range(len(pubmed)):
    # get the set of CNV genes for that study
    CNV_genes = get_human_CNV_genes_single_study(CNV_file, pubmed[i][0], 'all')
    print(pubmed[i][0], len(CNV_genes))
    # populate dict
    StudiesCNV[pubmed[i][0]] = CNV_genes
print('extracted CNV genes for each study')

# get synonym names for all genes {gene name : [list of synonyms]}
synonyms = get_synonyms('H_sapiens.gff3')
print('got gene synonyms', len(synonyms))    
# get the CDS sequences of the longest mRNAs for each gene {gene : sequence}
CDS_seq = extract_CDS_sequences('H_sapiens.gff3', 'H_sapiens_genome.txt', 'H_sapiens_valid_chromos.txt', keep_valid_chromos)
print('extracted CDS sequences', len(CDS_seq))  


# create a dict {study: {gene: CNV status}}    
CNV_status = {}
for study in StudiesCNV:
    # initialize dict
    CNV_status[study] = {}
    # loop over genes in CDS seq
    for gene in CDS_seq:
        # set boolean
        is_cnv = False
        # ask if gene in CNV genes
        if gene in StudiesCNV[study] or gene.upper() in StudiesCNV[study]:
            # gene is CNV, add gene and status to inner dict
            CNV_status[study][gene] = 'CNV'
        else:
            # ask if any of the gene synonyms are in CNV genes
            for name in synonyms[gene]:
                # check if gene in CNV
                if name in StudiesCNV[study] or name.upper() in StudiesCNV[study]:
                    # update boolean variable
                    is_cnv = True
            # check if gene in CNV
            if is_cnv == True:
                CNV_status[study][gene] = 'CNV'
            elif is_cnv == False:
                CNV_status[study][gene] = 'not_CNV'
print('sorted genes according to CNV status')    
    
    
# make a dictionary {study: {gene: [targets, seq_length, normalized_targets, CNV_status]}}
CNVTargets = {}
for study in CNV_status:
    # initialize inner dict
    CNVTargets[study] = {}
    # loop over genes with CNV status for given study
    for gene in CNV_status[study]:
        # populate with list of targets
        if gene in predicted_targets:
            CNVTargets[study][gene] = list(predicted_targets[gene])
            # add CNV status
            CNVTargets[study][gene].append(CNV_status[study][gene])
print('combined targets and CNV status')


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

    elif GraphType == 'bar':
        # Create a bar plot for proportions of replicates with CNV no diff on top of CNV lower
        ax.bar([0, 0.7, 1.4, 2.1], Data[3], width = 0.5, label = 'No difference', color= '#f7f7f7')
        # Create a bar plot for proportions of replicates with CNV greater on top of no diff
        ax.bar([0, 0.7, 1.4, 2.1], Data[0], width = 0.5, bottom = Data[3], label = 'CNV > non-CNV', color= '#e9a3c9')
        # Create a bar plot for proportions of replicates with CNV lower on top of CNV greater
        ax.bar([0, 0.7, 1.4, 2.1], Data[1], width = 0.5, bottom= Data[0], label = 'CNV < non-CNV', color = '#a1d76a')
    
      
    # write title   
    ax.set_title(Title, size = 8)
    
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
ax2 = CreateAx(2, 1, 2, Proportions, fig, figtitle, 0.45, Populations, [0.25, 0.95, 1.65, 2.35], 'bar')





## annotate Graph with significance level
#PvalTargetScan, PvalMiranda = [], []
#for species in species_names:
#    if CompTargetscan[species] >= 0.05:
#        PvalTargetScan.append('')
#    elif CompTargetscan[species] < 0.05 and CompTargetscan[species] >= 0.01:
#        PvalTargetScan.append('*')
#    elif CompTargetscan[species] < 0.01 and CompTargetscan[species] >= 0.001:
#        PvalTargetScan.append('**')
#    elif CompTargetscan[species] < 0.001:
#        PvalTargetScan.append('***')
#for species in species_names:
#    if CompMiranda[species] >= 0.05:
#        PvalMiranda.append('')
#    elif CompMiranda[species] < 0.05 and CompMiranda[species] >= 0.01:
#        PvalMiranda.append('*')
#    elif CompMiranda[species] < 0.01 and CompMiranda[species] >= 0.001:
#        PvalMiranda.append('**')
#    elif CompMiranda[species] < 0.001:
#        PvalMiranda.append('***')


## create list of Y and X positions to annotate figure with significance level
#if domain == 'CDS':
#    # make a list of Y positions
#    YposTargetscan = [0.41, 0.11, 0.16, 0.28, 0.14, 0.21]
#    YposMiranda = [0.32, 0.08, 0.11, 0.19, 0.11, 0.16]
#    Xpos = [0.2, 1.1, 2, 2.9, 3.8, 4.7]
#elif domain == '5UTR':
#    # make a list of Y positions
#    YposTargetscan = [0.41, 0.12, 0.17, 0.30, 0.16, 0.24]
#    YposMiranda = [0.32, 0.09, 0.12, 0.21, 0.12, 0.17]
#    Xpos = [0.2, 1.1, 2, 2.9, 3.8, 4.7]
#
#for i in range(len(PvalTargetScan)):
#    ax1.text(Xpos[i], YposTargetscan[i], PvalTargetScan[i], horizontalalignment = 'center',
#             verticalalignment = 'center', color = 'black', size = 8)
#for i in range(len(PvalMiranda)):
#    ax2.text(Xpos[i], YposMiranda[i], PvalMiranda[i], horizontalalignment = 'center',
#             verticalalignment = 'center', color = 'black', size = 8)

# make sure subplots do not overlap
plt.tight_layout()

## get outputfile
#if L == 15:
#    outputfile = 'PlotShort3UTR_15bp_' + domain + '_' + chromos + '_' + cnv_length + '_NormalizedSites' 
#elif L == 7:
#    outputfile = 'PlotShort3UTR_7bp_' + domain + '_' + chromos + '_' + cnv_length + '_NormalizedSites' 
#print(outputfile)

# save figure
fig.savefig('truc' + '.pdf', bbox_inches = 'tight')
       




















##############################

#
#
## use all chromos (including unplaced, unlocated, and MT) or only valid chromos 
## note: only CNVs on valid chromos are reported in DGV, so if all chromos are
## used it may introduce a bias by calling non CNV genes genes that cannot be detected
#
## keep genes on assembled nuclear chromosomes
#chromos = 'valid_chromos'
## consider all CNVs
#cnv_length = 'CNV_all_length'
#
## make a dictionary of species names : species code
#species_codes = {'H_sapiens': 'Hsa', 'P_troglodytes': 'Ptr', 'M_mulatta': 'Mml',
#                 'M_musculus': 'Mmu', 'B_taurus': 'Bta', 'G_gallus':'Gga'}
#
## make a dictionnary {species: {gene: [N_targets, Sequence_length, N_targets_normalized, CNV_status}}
## for each predictor, targetscan and miranda
#targetscan, miranda = {}, {}
#
#predictors = ['targetscan', 'miranda']
## loop over predictors
#for i in range(len(predictors)):
#    # loop over species names
#    for species in species_codes:
#        print(species)
#        
#        # get the seq input file
#        seq_input_file = species + '_' + domain + '_' + chromos + '_targetscan.txt'
#        print(seq_input_file)
#        
#        # get the predicted targets output file
#        predicted_targets = species + '_' + domain + '_' + chromos + '_predicted_sites_' + predictors[i] + '.txt'
#        print(predicted_targets)
#        
#        # parse the predictor outputfile to get a dict {gene: [targets, seq_length, normalized_targets]}
#        if i == 0:
#            targets = parse_targetscan_output(seq_input_file, predicted_targets, 'all')
#        elif i == 1:
#            targets = parse_miranda_output(seq_input_file, predicted_targets, 'all')
#        print('targets', len(targets))    
#        
#        # check if species is human
#        if species == 'H_sapiens':
#            # use DGV 2015 release 
#            CNV_file = 'H_sapiens_GRCh37_2015_CNV_all_length_valid_chromos.txt'
#        else:
#            CNV_file = species + '_' + cnv_length + '_' + chromos + '.txt' 
#        print(CNV_file)
#        
#        # get CNV gene status
#        CNV_status = sort_genes_CNV_status(CNV_file)
#        print('CNV status', len(CNV_status))
#        
#        # add CNV status
#        for gene in targets:
#            targets[gene].append(CNV_status[gene])
#            
#        # get file with UTR status (short, long)
#        UTR_file = species + '_3UTR_length_' + chromos + '.txt'
#        print(UTR_file)
#        # create a dict {gene : "short" (or "long")}
#        UTR_length = sort_genes_3UTR_length(UTR_file, L)
#        print('UTR length', len(UTR_length))
#        
#        # initialize inner dict
#        if i == 0:
#            targetscan[species] = {}
#            for gene in targets:
#                # record gene with short 3'UTR
#                if gene in UTR_length and UTR_length[gene] == 'short':
#                    targetscan[species][gene] = list(targets[gene])
#        elif i == 1:
#            miranda[species] = {}
#            for gene in targets:
#                # record gene with short 3'UTR
#                if gene in UTR_length and UTR_length[gene] == 'short':
#                    miranda[species][gene] = list(targets[gene])
#        
## check that the same number of species are recorded for miranda and targetscan
#assert len(targetscan) == len(miranda), 'different number of species depending on predictor'
## check that list values have the correct number of items
#for species in targetscan:
#    for gene in targetscan[species]:
#        assert len(targetscan[species][gene]) == 4, 'gene in targetscan does not have all required values'
#    for gene  in miranda[species]:
#        assert len(miranda[species][gene]) == 4, 'gene in miranda does not have all required values'
#
#
## create a dict of {species name : [[CNV], [non-CNV]]}
#SpeciesDataTargetscan, SpeciesDataMiranda = {}, {}
#
## loop over species names, get the number of normalized sites for CNV and non-CNV genes
#for species in targetscan:
#    # initialise list value
#    SpeciesDataTargetscan[species] = [[], []]
#    # populate inner lists with number of miRNA target sites per nucleotide
#    for gene in targetscan[species]:
#        if targetscan[species][gene][-1] == 'CNV':
#            SpeciesDataTargetscan[species][0].append(targetscan[species][gene][2])
#        elif targetscan[species][gene][-1] == 'not_CNV':
#            SpeciesDataTargetscan[species][1].append(targetscan[species][gene][2])
#for species in miranda:
#    # initialize list values
#    SpeciesDataMiranda[species] = [[], []]
#    # populate inner lists with number of mirna target sites per nucleotide
#    for gene in miranda[species]:
#        if miranda[species][gene][-1] == 'CNV':
#            SpeciesDataMiranda[species][0].append(miranda[species][gene][2])
#        elif miranda[species][gene][-1] == 'not_CNV':
#            SpeciesDataMiranda[species][1].append(miranda[species][gene][2])
#print('generated lists of target sites for CNV and non-CNV genes')
#
#
## perform stattistical tests between CNV and non-CNV genes
## create a dict to store results {species: P-value}
#CompTargetscan, CompMiranda = {}, {}
#for species in SpeciesDataTargetscan:
#    P = stats.ranksums(SpeciesDataTargetscan[species][0], SpeciesDataTargetscan[species][1])[1]
#    CompTargetscan[species] = P
#for species in SpeciesDataMiranda:
#    P = stats.ranksums(SpeciesDataMiranda[species][0], SpeciesDataMiranda[species][1])[1]    
#    CompMiranda[species] = P
#print('compared CNV and non-CNV genes')
#
## print P-values
#for species in CompTargetscan:
#    print('targetscan', species, CompTargetscan[species])
#for species in CompMiranda:
#    print('miranda', species, CompMiranda[species])
#
#
## make a list of data for each predictor
#AllDataTargetscan, AllDataMiranda = [], []
#
## make a list of species names to loop from
#species_names = ['H_sapiens', 'P_troglodytes', 'M_mulatta', 'M_musculus', 'B_taurus', 'G_gallus']
## loop over species in species names list and populate data lists, keeping the same order for targetscan and miranda
#for species in species_names:
#    # append list of target sites for CNV genes
#    AllDataTargetscan.append(SpeciesDataTargetscan[species][0])
#    AllDataMiranda.append(SpeciesDataMiranda[species][0])
#    # append list of target sites for non-CNV genes
#    AllDataTargetscan.append(SpeciesDataTargetscan[species][1])
#    AllDataMiranda.append(SpeciesDataMiranda[species][1])
#print('data consolidated in array')
#
#
## create figure
#fig = plt.figure(1, figsize = (8, 3))
#
## create list of labels and tick positions for the X axis
##xtickpos = [0.35, 1.25, 2.15, 3.05, 3.95, 4.85]
#xtickpos = [0.2, 1.1, 2, 2.9, 3.8, 4.7]
#Names = [species_codes[i] for i in species_names]
#print(Names)
#
## create a function to format the subplots
#def CreateAx(Columns, Rows, Position, Data, figure, Title, YMax, SpeciesNames, XScale):
#    '''
#    (int, int, int, list, figure_object, str, int, list, list)
#    Take the number of a column, and rows in the figure object and the position of
#    the ax in figure, a list of data, a title, a maximum value for the Y axis,
#    a list with species names and list of X axis tick positions and return an
#    ax instance in the figure
#    '''    
#    
#    
#    # create subplot in figure
#    # add a plot to figure (N row, N column, plot N)
#    ax = figure.add_subplot(Rows, Columns, Position)
#    # create a list of positions for the box plot    
#    BoxPositions = [0, 0.4, 0.9, 1.3, 1.8, 2.2, 2.7, 3.1, 3.6, 4, 4.5, 4.9]
#    # use a boxplot
#    bp = ax.boxplot(Data, showmeans = True, showfliers = False, widths = 0.3,
#                    positions = BoxPositions, patch_artist = True) 
#    
#    # color CNV and non-CNV boxes differently
#    i = 0    
#    # change box, whisker color to black
#    for box in bp['boxes']:
#        # change line color
#        box.set(color = 'black')
#        if i % 2 == 0:
#            # CNV data, color box in grey
#            box.set(facecolor = '#a6cee3')
#        else:
#            box.set(facecolor = '#b2df8a')
#        i += 1
#    # change whisker color to black
#    for wk in bp['whiskers']:
#        wk.set(color = 'black', linestyle = '-')
#    # change color of the caps
#    for cap in bp['caps']:
#        cap.set(color = 'black')
#    # change the color and line width of the medians
#    for median in bp['medians']:
#        median.set(color = 'black')
#    # change the mean marker and marker
#    for mean in bp['means']:
#        mean.set(marker = 'o', markeredgecolor = 'black', markerfacecolor = 'black', markersize = 3)
#    
#    
#    # write title   
#    ax.set_title(Title, size = 8)
#    
#    # set font for all text in figure
#    FigFont = {'fontname':'Arial'}   
#    
#    # write label for y axis
#    ax.set_ylabel('Normalized number of miRNA\nsites per gene', color = 'black',  size = 8, ha = 'center', **FigFont)
#
#    # write label for x axis
#    plt.xticks(XScale, SpeciesNames, ha = 'center', fontsize = 8, **FigFont)
#
#    # add a range for the Y axis
#    plt.ylim([0, YMax])
#    
#    plt.xlim([-0.25, 5.15])
#    
#    # do not show lines around figure  
#    ax.spines["top"].set_visible(False)    
#    ax.spines["bottom"].set_visible(True)    
#    ax.spines["right"].set_visible(False)    
#    ax.spines["left"].set_visible(False)  
#    # offset the spines
#    for spine in ax.spines.values():
#        spine.set_position(('outward', 5))
#    
#    # add a light grey horizontal grid to the plot, semi-transparent, 
#    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5, linewidth = 0.5)  
#    # hide these grids behind plot objects
#    ax.set_axisbelow(True)
#
#    # do not show ticks
#    plt.tick_params(
#        axis='both',       # changes apply to the x-axis and y-axis (other option : x, y)
#        which='both',      # both major and minor ticks are affected
#        bottom='on',      # ticks along the bottom edge are off
#        top='off',         # ticks along the top edge are off
#        right = 'off',
#        left = 'off',          
#        labelbottom='on', # labels along the bottom edge are on
#        colors = 'black',
#        labelsize = 8,
#        direction = 'out') # ticks are outside the frame when bottom = 'on'  
#      
#    # Set the tick labels font name
#    for label in ax.get_yticklabels():
#        label.set_fontname('Arial')
#    
#    # create a margin around the x axis
#    plt.margins(0.05)
#    
#    return ax      
#
#
## plot data for targetscan
#ax1 = CreateAx(2, 1, 1, AllDataTargetscan, fig, 'TargetScan', 0.45, Names, xtickpos)
#ax2 = CreateAx(2, 1, 2, AllDataMiranda, fig, 'miRanda', 0.45, Names, xtickpos)
#
## annotate Graph with significance level
#PvalTargetScan, PvalMiranda = [], []
#for species in species_names:
#    if CompTargetscan[species] >= 0.05:
#        PvalTargetScan.append('')
#    elif CompTargetscan[species] < 0.05 and CompTargetscan[species] >= 0.01:
#        PvalTargetScan.append('*')
#    elif CompTargetscan[species] < 0.01 and CompTargetscan[species] >= 0.001:
#        PvalTargetScan.append('**')
#    elif CompTargetscan[species] < 0.001:
#        PvalTargetScan.append('***')
#for species in species_names:
#    if CompMiranda[species] >= 0.05:
#        PvalMiranda.append('')
#    elif CompMiranda[species] < 0.05 and CompMiranda[species] >= 0.01:
#        PvalMiranda.append('*')
#    elif CompMiranda[species] < 0.01 and CompMiranda[species] >= 0.001:
#        PvalMiranda.append('**')
#    elif CompMiranda[species] < 0.001:
#        PvalMiranda.append('***')
#
#
## create list of Y and X positions to annotate figure with significance level
#if domain == 'CDS':
#    # make a list of Y positions
#    YposTargetscan = [0.41, 0.11, 0.16, 0.28, 0.14, 0.21]
#    YposMiranda = [0.32, 0.08, 0.11, 0.19, 0.11, 0.16]
#    Xpos = [0.2, 1.1, 2, 2.9, 3.8, 4.7]
#elif domain == '5UTR':
#    # make a list of Y positions
#    YposTargetscan = [0.41, 0.12, 0.17, 0.30, 0.16, 0.24]
#    YposMiranda = [0.32, 0.09, 0.12, 0.21, 0.12, 0.17]
#    Xpos = [0.2, 1.1, 2, 2.9, 3.8, 4.7]
#
#for i in range(len(PvalTargetScan)):
#    ax1.text(Xpos[i], YposTargetscan[i], PvalTargetScan[i], horizontalalignment = 'center',
#             verticalalignment = 'center', color = 'black', size = 8)
#for i in range(len(PvalMiranda)):
#    ax2.text(Xpos[i], YposMiranda[i], PvalMiranda[i], horizontalalignment = 'center',
#             verticalalignment = 'center', color = 'black', size = 8)
#
## make sure subplots do not overlap
#plt.tight_layout()
#
## get outputfile
#if L == 15:
#    outputfile = 'PlotShort3UTR_15bp_' + domain + '_' + chromos + '_' + cnv_length + '_NormalizedSites' 
#elif L == 7:
#    outputfile = 'PlotShort3UTR_7bp_' + domain + '_' + chromos + '_' + cnv_length + '_NormalizedSites' 
#print(outputfile)
#
## save figure
#fig.savefig(outputfile + '.eps', bbox_inches = 'tight')
#       
