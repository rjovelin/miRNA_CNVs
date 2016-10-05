# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 16:21:29 2016

@author: RJovelin
"""


# use this script to compare target sites for short 3' UTRs genes using 
# Zarrei's CNV map of the human genome

# This script keeps genes on assembled nuclear chromosomes and considers all CNVs
# and uses the stringent CNV map of the human genome

# usage PlotTargetShort3UTRGoldCNVMap [options]
# [pdf/ai/png]: save as eps if no argument is provided or as the format indicated by the argument

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
import copy
# import custom modules
from CNV_miRNAs import *


if len(sys.argv) == 1:
    extension = '.eps'
elif len(sys.argv) == 2:
    extension = sys.argv[1]
    assert extension in ['pdf', 'ai', 'png']
    extension = '.' + extension


# use this function to get the targets per gene per domain
# for genes with short 3'UTRs {domain: {gene: [targets, seq_length, normalized_targets, CNV status]}}
def GetTargetsPerGeneDomain(predictor, L, GeneCNV, UTR_file = 'H_sapiens_3UTR_length_valid_chromos.txt'):
    '''
    (str, int, dict) -> dict
    Take the name of the predictor algorithm, the length threshold for short 3'UTRs,
    a dictionary with CNV gene status and return a dictionary with target sites
    per nucleotide and per gene in 5'UTRs and CDS of genes with short 3'UTRs    
    '''
    
    # create a dict with targets per region
    # {domain: {gene: [targets, seq_length, normalized_targets, CNV status]}}
    TargetSites = {}
    
    # create a dict {gene : "short" (or "long")}
    UTR_length = sort_genes_3UTR_length(UTR_file, L)
    print('UTR length', len(UTR_length))
    
    # loop over 5'UTR and CDS regions
    for domain in ['5UTR', 'CDS']:
        # get the seq input file
        seq_input_file = 'H_sapiens_' + domain + '_valid_chromos_targetscan.txt'
        print(seq_input_file)
        # get the predicted targets output file
        predicted_targets = 'H_sapiens_' + domain + '_valid_chromos_predicted_sites_' + predictor + '.txt'
        print(predicted_targets)
        # parse the predictor outputfile to get a dict {gene: [targets, seq_length, normalized_targets]}
        if predictor == 'targetscan':
            targets = parse_targetscan_output(seq_input_file, predicted_targets, 'all')
        elif predictor == 'miranda':
            targets = parse_miranda_output(seq_input_file, predicted_targets, 'all')
        print('targets', len(targets))    
    
        # add CNV status
        for gene in targets:
            if gene in GeneCNV:
                targets[gene].append(GeneCNV[gene])
        
        # initialize inner dict
        TargetSites[domain] = {}
        for gene in targets:
            # record gene with short 3'UTR
            if gene in UTR_length and UTR_length[gene] == 'short':
                TargetSites[domain][gene] = copy.deepcopy(targets[gene])
    # remove genes if genes do not have CNV status
    for domain in TargetSites:
        to_remove, to_keep = [], []
        for gene in TargetSites[domain]:
            if len(TargetSites[domain][gene]) != 4:
                assert len(TargetSites[domain][gene]) == 3, 'gene without CNV status should have correct target information'
                assert type(TargetSites[domain][gene][-1]) != str , 'last item in the list should be a number if CNV status is missing'
                to_remove.append(gene)
            else:
                to_keep.append(gene)
        # remove genes without a CNV status
        print('remove {0} genes from {1}'.format(len(to_remove), domain))    
        for gene in to_remove:
            del TargetSites[domain][gene]
        assert len(TargetSites[domain]) == len(to_keep), 'numbers of genes with CNV status and targetscan targets should match'

    return TargetSites


# use this function to group targets in lists
# [[5UTR targets cnv], [5UTR targets noncnv], [CDS targets cnv], [CDS targets non-cnv]]
def TargetsToArray(TargetSites):
    '''
    (dict) -> dict
    Take the dictionary with targets per gene and domain and return a list with
    with targets in domains for CNV and non-CNV genes with short 3'UTRs
    '''
    # create a list of lists of targets for CNV and non-CNV genes for each region
    # [[5UTR targets cnv], [5UTR targets noncnv], [CDS targets cnv], [CDS targets non-cnv]]
    AllData = [[], [], [], []]
    for domain in TargetSites:
        for gene in TargetSites[domain]:
            if domain == '5UTR' and TargetSites[domain][gene][-1] == 'CNV':
                i = 0
            elif domain == '5UTR' and TargetSites[domain][gene][-1] == 'not_CNV':
                i = 1
            elif domain == 'CDS' and TargetSites[domain][gene][-1] == 'CNV':
                i = 2
            elif domain == 'CDS' and TargetSites[domain][gene][-1] == 'not_CNV':
                i = 3
            AllData[i].append(TargetSites[domain][gene][2])        
    return AllData


# use this function to perform stattistical tests between CNV and non-CNV genes
# create a list to store significance level [list of significance levels]
def GetSignificanceLevel(AllData):
    '''
    (list) -> list
    Take the list of targets for domains and CNV status and return a list
    with significance level
    '''
    
    Significance = []
    # compare CNV and non-CNV genes
    for i in range(0, len(AllData), 2):
        Pval = stats.ranksums(AllData[i], AllData[i+1])[1]
        if Pval >= 0.05:
            Significance.append('')
        elif Pval < 0.05 and Pval >= 0.01:
            Significance.append('*')
        elif Pval < 0.01 and Pval >= 0.001:
            Significance.append('**')
        elif Pval < 0.001:
            Significance.append('***')
    return Significance


# create a function to format the subplots
def CreateAx(Columns, Rows, Position, Data, figure, Title, Ymax, YAxisLine):
    '''
    (int, int, int, list, figure_object, str, str, float, bool)
    Take the number of a column, and rows in the figure object and the position of
    the ax in figure, the predictor algorithm used, a list of data, a title,
    the maximum value of the Y axis, a boolean to specify if Y axis should be
    specified and return an ax instance in the figure
    '''    
    
    # create subplot in figure
    # add a plot to figure (N row, N column, plot N)
    ax = figure.add_subplot(Rows, Columns, Position)
    # create a list of positions for the box plot    
    BoxPositions = [0, 0.4, 0.9, 1.3]

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
            # CNV data
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

    # write title   
    ax.set_title(Title, size = 8)
    # set font for all text in figure
    FigFont = {'fontname':'Arial'}   
    # write label for x and y axis
    ax.set_ylabel('miRNA sites / nt', color = 'black',  size = 8, ha = 'center', **FigFont)
    if YAxisLine == True:
        # write label for x axis
        plt.xticks([0.2, 1.1], ['5\'UTR', 'CDS'], ha = 'center', fontsize = 8, **FigFont)
    
    # add a range for the Y axis
    plt.ylim([0, Ymax])
    
    # do not show lines around figure  
    ax.spines["top"].set_visible(False)    
    ax.spines["right"].set_visible(False)    
    ax.spines["left"].set_visible(False)  
    if YAxisLine == False:
        ax.spines["bottom"].set_visible(False)    
    elif YAxisLine == True:
        ax.spines["bottom"].set_visible(True)    
        # offset the spines
        for spine in ax.spines.values():
            spine.set_position(('outward', 5))
    
    # add a light grey horizontal grid to the plot, semi-transparent, 
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5, linewidth = 0.5)  
    # hide these grids behind plot objects
    ax.set_axisbelow(True)

    if YAxisLine == True:
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
    elif YAxisLine == False:
        # do not show ticks
        plt.tick_params(
            axis='both',       # changes apply to the x-axis and y-axis (other option : x, y)
            which='both',      # both major and minor ticks are affected
            bottom='off',      # ticks along the bottom edge are off
            top='off',         # ticks along the top edge are off
            right = 'off',
            left = 'off',          
            labelbottom='off', # labels along the bottom edge are on
            colors = 'black',
            labelsize = 8,
            direction = 'out') # ticks are outside the frame when bottom = 'on'  
    
    # Set the tick labels font name
    for label in ax.get_yticklabels():
        label.set_fontname('Arial')
    # create a margin around the x axis
    plt.margins(0.05)
    
    return ax      




# make a list of fasta files
files = [i for i in os.listdir('./GRCH37_genome') if i[-3:] == '.fa']
print(len(files))

# make a dictionary of scaffold: chromosome
chromos = {}
# loop over fasta files
for filename in files:
    infile = open('./GRCH37_genome/' + filename)
    # look for sequence headers
    for line in infile:
        if line.startswith('>'):
            line = line.rstrip().split('|')
            scaffold = line[3]
            LG = line[-1]
            LG = LG[:LG.index(',')]
            LG = LG.replace('Homo sapiens', '')
            LG = LG.replace(' ', "")
            LG = LG.replace('omosome', '')
            assert scaffold not in chromos, 'scaffold is already recorded'            
            chromos[scaffold] = LG
    infile.close()
print('done matching chromosome names', len(chromos))
                
# get CNVR coordinates {CNVR: [chromo, start, end]}
# use stringent map
CNVFile = './GRCH37_genome/Stringent.Gain+Loss.hg19.2015-02-03.txt'

CNVCoord = {}
infile = open(CNVFile)   
# skip header, read file
infile.readline()
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split()
        # get chromo, start, end, state, region
        chromo, start, end, state, CNVR = line[0], int(line[1]) -1, int(line[2]), line[3], line[4]
        # check that state is CNV
        if state == 'CNV':
            CNVCoord[CNVR] = [chromo, start, end]
infile.close()
print('got CNVR coordinates', len(CNVCoord))


# record CNVR per chromo
CNVRChromo = {}
for cnv in CNVCoord:
    chromo = CNVCoord[cnv][0]
    if chromo in CNVRChromo:
        CNVRChromo[chromo].append(cnv)
    else:
        CNVRChromo[chromo] = [cnv]
print('recorded CNVR for each chromosome')

# get the gene coordinates, keeping the longest mRNA per gene
GFF_file = './GRCH37_genome/ref_GRCh37.p5_top_level.gff3'       
       
# match gene ID to gene name  {gene ID : gene name}
GeneIDToGeneName = {}
# open file for reading
infile = open(GFF_file)
# loop over file
for line in infile:
    # get line with mRNA
    if 'gene' in line:
        line = line.rstrip().split('\t')
        if line[2] == 'gene':
            # get chromo
            chromo = line[0]
            # parse descriptor string
            # separate on ';' because line doesn't always have the same structure
            description = line[-1].split(';')
            # loop over strings in list, find and extract GeneID
            for i in range(len(description)):
                if 'GeneID:' in description[i]:
                    ID = description[i]
                if 'Name=' in description[i]:
                    name = description[i]
            # parse ID and name
            if ',' in ID:
                # check if comma happens before or after geneID
                if ID.count(',') == 1 and ID.index(',') < ID.index('GeneID:'):
                    ID = ID[ID.index('GeneID:') + 7:]
                else:
                    ID = ID[ID.index('GeneID:') + 7: ID.index(',', ID.index('GeneID:'))]
            else:
                ID = ID[ID.index('GeneID:') + 7:]
            if ',' in name:
                name = name[name.index('Name=') + 5: name.index(',', name.index('Name='))]
            else:
                name = name[name.index('Name=') + 5:]
            GeneIDToGeneName[ID] = name
# close file
infile.close()
print('matched gene ID to gene name', len(GeneIDToGeneName))

# match RNA ID to gene ID {RNA ID : gene ID} 
mRNAToGene = {}
# open file for reading
infile = open(GFF_file, 'r')
# loop over file
for line in infile:
    # find lines with mRNA
    if 'mRNA' in line:
        line = line.rstrip().split('\t')
        if line[2] == 'mRNA':
            # extract rna id
            rna_id = line[-1][line[-1].index('ID=') + 3: line[-1].index(';')]
            # parse description line and extract gene ID
            description = line[-1].split(';')
            for i in range(len(description)):
                # find and extract gene ID
                if 'GeneID:' in description[i]:
                    gene = description[i]
            # further parse transcript ID and gene ID
            if ',' in gene:
                # check if comma happens before or after geneID
                if gene.count(',') == 1 and gene.index(',') < gene.index('GeneID:'):
                    gene = gene[gene.index('GeneID:') + 7:]
                else:
                    gene = gene[gene.index('GeneID:') + 7: gene.index(',', gene.index('GeneID:'))]                    
            else:
                gene = gene[gene.index('GeneID:') + 7: ]                    
            # populate dict 
            mRNAToGene[rna_id] = gene
# close file
infile.close()
print('matched mRNA ID to gene ID', len(mRNAToGene))

# make a dict with gene ID as key and a list of corresponding mRNA ID {gene ID: [rna ID]}
GeneTomRNA = {}
for rna in mRNAToGene:
    if mRNAToGene[rna] in GeneTomRNA:
        GeneTomRNA[mRNAToGene[rna]].append(rna)
    else:
        GeneTomRNA[mRNAToGene[rna]] = [rna]
print('matched Gene ID with all their mRNA ID', len(GeneTomRNA))

# get the coordinates of all mRNAs
# create a dict {rna_id: [chromo, start, end , orientation]}
mRNACoord = {}
# open file for reading
infile = open(GFF_file, 'r')
# loop over file
for line in infile:
    # ignore lines that do not contain mRNA
    if 'mRNA' in line:
        line = line.rstrip().split('\t')
        # check that line correspond to a a mRNA
        if line[2] == 'mRNA':
            # extract RNA ID
            rna_id = line[-1][line[-1].index('ID=') + 3: line[-1].index(';')]
            # get chromo
            LG = line[0]
            # get chromo name (eg chr1)
            if LG in chromos:
                chromo = chromos[LG]
                # get orientation
                orientation = line[6]
                # get start, end positions 0-based
                start = int(line[3]) -1
                end = int(line[4])
                mRNACoord[rna_id] = [chromo, start, end, orientation]  
# close file
infile.close()
print('extracted mRNA coordinates', len(mRNACoord))    

# remove rna with no coord
for gene in GeneTomRNA:
    to_remove = []
    for rna in GeneTomRNA[gene]:
        if rna not in mRNACoord:
            to_remove.append(rna)
    print('rna to remove: {0}'.format(len(to_remove)), end = '\r')
    for rna in to_remove:
        GeneTomRNA[gene].remove(rna)
print('removed rna without coordinates')
# remove genes without rnas
to_remove = [gene for gene in GeneTomRNA if len(GeneTomRNA[gene]) == 0]
if len(to_remove) != 0:
    for gene in to_remove:
        del GeneTomRNA[gene]
    print('removed genes without rnas', len(GeneTomRNA))

# remove the few genes that are mapped to different chromosomes
# match the chromos of each mRNA to their parent gene {gene: {set of chromos}}
checkgene = {}
for gene in GeneTomRNA:
    # get the chromo of each rna
    for rna in GeneTomRNA[gene]:
        LG = mRNACoord[rna][0]
        if gene not in checkgene:
            checkgene[gene] = set()
        checkgene[gene].add(LG)
# remove genes mapped to more than 1 chromo
to_remove = []
for gene in GeneTomRNA:
    if len(checkgene[gene]) > 1:
        to_remove.append(gene)            
for gene in to_remove:
    del GeneTomRNA[gene]
print('removed genes mapped to more than 1 chromosome', len(GeneTomRNA))

# record gene per chromosome
GeneChromo = {}
for gene in GeneTomRNA:
    # get the chromo of the first rna
    chromo = mRNACoord[GeneTomRNA[gene][0]][0]
    # loop over the rna of that gene
    for rna in GeneTomRNA[gene]:
        LG = mRNACoord[rna][0]
        assert chromo == LG, 'mRNAs of the same gene should be on the same chromosome'
    if chromo in GeneChromo:
        GeneChromo[chromo].append(gene)
    else:
        GeneChromo[chromo] = [gene]
print('recorded genes for each chromosome')
       
# count the number of cnv and non-cnv mRNAs
a, b, c = 0, 0, len(GeneTomRNA)
# record the CNV status of all genes {gene name: CNV status}
GeneCNV = {}
# find all genes affected by CNVR
for chromo in GeneChromo:
    # loop through gene on that chromo
    for gene in GeneChromo[chromo]:
        # get gene name
        name = GeneIDToGeneName[gene]
        # update counter
        c -= 1
        # loop over the gene's mRNAs
        for rna in GeneTomRNA[gene]:
            # set boolean
            FoundCNV = False
            # get mRNA coord
            rna_chromo, rna_start, rna_end = mRNACoord[rna][0], mRNACoord[rna][1], mRNACoord[rna][2]
            # loop through CNVR on that chromo
            for CNVR in CNVRChromo[chromo]:
                # get chromo, start and end
                cnv_chromo, cnv_start, cnv_end  = CNVCoord[CNVR][0], CNVCoord[CNVR][1], CNVCoord[CNVR][2]
                assert cnv_chromo == rna_chromo, 'chromos for cnv and rna should match'
                # check positions to see if rna is affected by CNVR                
                if rna_start < cnv_end and (rna_end > cnv_start or rna_end > cnv_end):
                    FoundCNV = True
                    # exit loop, no need to check the other CNVR
                    break
            # check if the mrna overlaps with a CNV region
            if FoundCNV == True:
                # rna is found in CNVR, update gene status
                GeneCNV[name] = 'CNV'
                # update counter and exit loop, no need to check other mRNAs
                a += 1
                break
        # check if any rna has been found in CNVR
        if FoundCNV == False:
            assert name not in GeneCNV, 'CNV status for that gene should not have been already recorded'
            # update counter and CNV status after all the mRNAs have beeb checked
            GeneCNV[name] = 'not_CNV'
            b += 1
        print('chromo: {0}, cnv: {1}, non-cnv: {2}, remaining: {3}'.format(chromo, a, b, c), sep = '\t', end = '\r')   

print('\n')
print('chromo: {0}, cnv: {1}, non-cnv: {2}, remaining: {3}'.format(chromo, a, b, c), sep = '\t', end = '\n')   
             

# make sets of cnv and noncnv genes, including synonyms
# add all synonymous genes
cnv, noncnv = set(), set()
for gene in GeneCNV:
    if GeneCNV[gene] == 'CNV':
        cnv.add(gene)
    elif GeneCNV[gene] == 'not_CNV':
        noncnv.add(gene)
assert len(cnv.intersection(noncnv)) == 0, 'genes cannot be both CNV and non-CNV'


# get the targets per gene and domain
TargetScanTargets7 = GetTargetsPerGeneDomain('targetscan', 7, GeneCNV, UTR_file = 'H_sapiens_3UTR_length_valid_chromos.txt')
TargetScanTargets15 = GetTargetsPerGeneDomain('targetscan', 15, GeneCNV, UTR_file = 'H_sapiens_3UTR_length_valid_chromos.txt')
MirandaTargets7 = GetTargetsPerGeneDomain('miranda', 7, GeneCNV, UTR_file = 'H_sapiens_3UTR_length_valid_chromos.txt')
MirandaTargets15 = GetTargetsPerGeneDomain('miranda', 15, GeneCNV, UTR_file = 'H_sapiens_3UTR_length_valid_chromos.txt')
print('got targets per gene and domain')

# check that all domains have been recorded for miranda and targetscan
assert len(TargetScanTargets7) == len(TargetScanTargets15) == len(MirandaTargets15) == len(MirandaTargets7) == 2, 'different number of domains depending on predictor'


# group targets in list for each domain and CNV status
# [[5UTR targets cnv], [5UTR targets noncnv], [CDS targets cnv], [CDS targets non-cnv]]
AllDataTargetScan7 = TargetsToArray(TargetScanTargets7)
AllDataTargetScan15 = TargetsToArray(TargetScanTargets15)
AllDataMiranda7 = TargetsToArray(MirandaTargets7)
AllDataMiranda15 = TargetsToArray(MirandaTargets15)
print('data consolidated in array')


# Perform stattistical tests between CNV and non-CNV genes [list of significance levels]
SignificanceTargetscan7 = GetSignificanceLevel(AllDataTargetScan7)
SignificanceTargetscan15 = GetSignificanceLevel(AllDataTargetScan15)
SignificanceMiranda7 = GetSignificanceLevel(AllDataMiranda7)
SignificanceMiranda15 = GetSignificanceLevel(AllDataMiranda15)
print('compared CNV and non-CNV genes')


# plot the number of targets for CNV and non-CNV genes for each region

# create list of labels and tick positions for the X axis
#xtickpos = [0.2, 1.1, 2, 2.9]


# create figure
fig = plt.figure(1, figsize = (3.4, 3))
# plot data 
ax1 = CreateAx(2, 2, 1, AllDataTargetScan7, fig, 'Targetscan', 0.45, False)
ax2 = CreateAx(2, 2, 2, AllDataMiranda7, fig, 'miRanda', 0.45,  False)
ax3 = CreateAx(2, 2, 3, AllDataTargetScan15, fig, 'Targetscan', 0.45, True)
ax4 = CreateAx(2, 2, 4, AllDataMiranda15, fig, 'miRanda', 0.45, True)


# annotate figure with significance level
# create list of Y and X positions to annotate figure with significance level
Xpos = [0.2, 1.1]
Ypos_targetscan = [0.42, 0.40]
Ypos_miranda = [0.35, 0.32]

for i in range(len(SignificanceTargetscan7)):
    ax1.text(Xpos[i], Ypos_targetscan[i], SignificanceTargetscan7[i], horizontalalignment = 'center', verticalalignment = 'center', color = 'black', size = 8)
    ax2.text(Xpos[i], Ypos_miranda[i], SignificanceMiranda7[i], horizontalalignment = 'center', verticalalignment = 'center', color = 'black', size = 8)
    ax3.text(Xpos[i], Ypos_targetscan[i], SignificanceTargetscan15[i], horizontalalignment = 'center', verticalalignment = 'center', color = 'black', size = 8)
    ax4.text(Xpos[i], Ypos_miranda[i], SignificanceMiranda15[i], horizontalalignment = 'center', verticalalignment = 'center', color = 'black', size = 8)


# add legend relative to ax1 using ax1 coordinates
C = mpatches.Patch(facecolor = '#a6cee3', edgecolor = 'black', linewidth = 1, label= 'CNV')
N = mpatches.Patch(facecolor = '#b2df8a', edgecolor = 'black', linewidth = 1, label= 'non-CNV')
ax1.legend(handles = [C, N], loc = (-0.25, 1.2), fontsize = 8, frameon = False, ncol = 2)

# make sure subplots do not overlap
plt.tight_layout()

# add subplot label
ax1.text(-1.4, 0.51, 'A', horizontalalignment = 'center',
         verticalalignment = 'center', color = 'black', size = 8)
ax3.text(-1.4, 0.51, 'B', horizontalalignment = 'center',
         verticalalignment = 'center', color = 'black', size = 8)   

# build outputfile with arguments
outputfile = 'PlotTargetShort3UTRGoldMapStringent_CNVAllLength_ValidChromos'
print(outputfile)

# save figure
fig.savefig(outputfile + extension, bbox_inches = 'tight')    

