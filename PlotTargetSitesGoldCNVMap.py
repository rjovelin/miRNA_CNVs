# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 14:21:28 2016

@author: RJovelin
"""

# use this script to make a file of chromosome names matched to scaffold/contig names

# usage PlotTargetSitesGoldCNVMap [options]
# - [stringent/inclusive]: filter level to define CNVR


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



# get option to use stingent or inclusive CNV definitions
CNVFilter = sys.argv[1]
assert CNVFilter in ['stringent', 'inclusive'], 'should use appropriate option'

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
if CNVFilter == 'stringent':
    CNVFile = './GRCH37_genome/Stringent.Gain+Loss.hg19.2015-02-03.txt'
elif CNVFilter == 'inclusive':
    CNVFile = './GRCH37_genome/Inclusive.Gain+Loss.hg19.2015-02-03.txt'

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


# make a dictionary with domain as key and a dictionary of gene: normalized targets as value
# {domain: {gene: [targets, seq_length, normalized_targets]}}
regions = ['3UTR', '5UTR', 'CDS']
targetscan, miranda = {}, {} 
for domain in regions:
    # get the seq input file
    seq_input_file = 'H_sapiens_' + domain + '_valid_chromos_targetscan.txt'
    # get the predicted targets output file {gene: [targets, seq_length, normalized_targets]}
    predicted_targets = 'H_sapiens_' + domain + '_valid_chromos_predicted_sites_targetscan.txt'
    targets = parse_targetscan_output(seq_input_file, predicted_targets, 'all')
    targetscan[domain] = dict(targets)
print('got targetscan targets for each domain of each gene')
for domain in regions:
    # get the seq input file
    seq_input_file = 'H_sapiens_' + domain + '_valid_chromos_targetscan.txt'
    # get the predicted targets output file {gene: [targets, seq_length, normalized_targets]}
    predicted_targets = 'H_sapiens_' + domain + '_valid_chromos_predicted_sites_miranda.txt'
    targets = parse_miranda_output(seq_input_file, predicted_targets, 'all')
    miranda[domain] = dict(targets)
print('got miranda targets for each domain of each gene')


# remove genes if genes do not have CNV status
for region in targetscan:
    to_remove, to_keep = [], []
    for gene in targetscan[region]:
        if gene in cnv or gene in noncnv:
            to_keep.append(gene)
        elif gene not in cnv and gene not in noncnv:
            to_remove.append(gene)
    # remove genes without a CNV status
    print('remove {0} genes from {1}'.format(len(to_remove), region))    
    for gene in to_remove:
        del targetscan[region][gene]
    assert len(targetscan[region]) == len(to_keep), 'numbers of genes with CNV status and targetscan targets should match'
for region in miranda:
    to_remove, to_keep = [], []
    for gene in miranda[region]:
        if gene in cnv or gene in noncnv:
            to_keep.append(gene)
        elif gene not in cnv and gene not in noncnv:
            to_remove.append(gene)
    # remove genes without a CNV status
    print('remove {0} genes from {1}'.format(len(to_remove), region))  
    for gene in to_remove:
        del miranda[region][gene]
    assert len(miranda[region]) == len(to_keep), 'numbers of genes with CNV status and miranda targets should match'
print('removed genes without CNV status')

# add CNV status of genes with target prediction
for region in targetscan:
    for gene in targetscan[region]:
        if gene in cnv:
            targetscan[region][gene].append('CNV')
        elif gene in noncnv:
            targetscan[region][gene].append('not_CNV')
for region in miranda:
    for gene in miranda[region]:
        if gene in cnv:
            miranda[region][gene].append('CNV')
        elif gene in noncnv:
            miranda[region][gene].append('not_CNV')
print('added CNV status to each gene domain')    
    
 
# plot the number of targets for CNV and non-CNV genes for each region

# create a dict with predictor as key and a list of targets for CNV and non-CNV genes for each region
# {predictor: [[3UTR targets cnv], [3UTR targets non-cnv], [5UTR targets cnv], [5UTR targets noncnv], [CDS targets cnv], [CDS targets non-cnv]]
AllData = {}
AllData['targetscan'], AllData['miranda'] = [[], [], [], [], [], []], [[], [], [], [], [], []]
for region in targetscan:
    for gene in targetscan[region]:
        if region == '3UTR' and targetscan[region][gene][-1] == 'CNV':
            i = 0
        elif region == '3UTR' and targetscan[region][gene][-1] == 'not_CNV':
            i = 1
        elif region == '5UTR' and targetscan[region][gene][-1] == 'CNV':
            i = 2
        elif region == '5UTR' and targetscan[region][gene][-1] == 'not_CNV':
            i = 3
        elif region == 'CDS' and targetscan[region][gene][-1] == 'CNV':
            i = 4
        elif region == 'CDS' and targetscan[region][gene][-1] == 'not_CNV':
            i = 5
        AllData['targetscan'][i].append(targetscan[region][gene][2])        
for region in miranda:
    for gene in miranda[region]:
        if region == '3UTR' and miranda[region][gene][-1] == 'CNV':
            i = 0
        elif region == '3UTR' and miranda[region][gene][-1] == 'not_CNV':
            i = 1
        elif region == '5UTR' and miranda[region][gene][-1] == 'CNV':
            i = 2
        elif region == '5UTR' and miranda[region][gene][-1] == 'not_CNV':
            i = 3
        elif region == 'CDS' and miranda[region][gene][-1] == 'CNV':
            i = 4
        elif region == 'CDS' and miranda[region][gene][-1] == 'not_CNV':
            i = 5
        AllData['miranda'][i].append(miranda[region][gene][2])    
print('data consolidated in array')

for predictor in AllData:
    print(predictor, end = ' ')
    print(' '.join(list(map(lambda x: str(len(x)), AllData[predictor]))))     

# perform stattistical tests between CNV and non-CNV genes
# create dicts to store results {predictor: [list of significance levels]}
Significance = {}
for predictor in AllData:
    Significance[predictor] = []
    # compare CNV and non-CNV genes
    for i in range(0, len(AllData[predictor]), 2):
        Pval = stats.ranksums(AllData[predictor][i], AllData[predictor][i+1])[1]
        if Pval >= 0.05:
            Significance[predictor].append('')
        elif Pval < 0.05 and Pval >= 0.01:
            Significance[predictor].append('*')
        elif Pval < 0.01 and Pval >= 0.001:
            Significance[predictor].append('**')
        elif Pval < 0.001:
            Significance[predictor].append('***')
print('compared CNV and non-CNV genes')


# create list of labels and tick positions for the X axis
#xtickpos = [0.2, 1.1, 2, 2.9, 3.8, 4.7]

# create a function to format the subplots
def CreateAx(Columns, Rows, Position, predictor, Data, figure, Title, YAxisLine):
    '''
    (int, int, int, str, list, figure_object, str, str, bool)
    Take the number of a column, and rows in the figure object and the position of
    the ax in figure, the predictor algorithm used, a list of data, a title,
    a boolean to specify if Y axis should be specified and return an
    ax instance in the figure
    '''    
    
    # create subplot in figure
    # add a plot to figure (N row, N column, plot N)
    ax = figure.add_subplot(Rows, Columns, Position)
    # create a list of positions for the box plot    
    BoxPositions = [0, 0.4, 0.9, 1.3, 1.8, 2.2]

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
        plt.xticks([0.2, 1.1, 2], ['3\'UTR', '5\'UTR', 'CDS'], ha = 'center', fontsize = 8, **FigFont)
    # add a range for the Y axis
    if predictor == 'targetscan':
        plt.ylim([0.10, 0.46])
    elif predictor == 'miranda':
        plt.ylim([0, 0.36])
    plt.xlim([-0.25, 2.45])

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

# create figure
fig = plt.figure(1, figsize = (2.3, 3.1))
# plot data 
ax1 = CreateAx(1, 2, 1, 'targetscan', AllData['targetscan'], fig, 'Targetscan', False)
ax2 = CreateAx(1, 2, 2, 'miranda', AllData['miranda'], fig, 'miRanda', True)


# annotate figure with significance level
# create list of Y and X positions to annotate figure with significance level
Xpos = [0.2, 1.1, 2]
Ypos_targetscan = [0.42, 0.40, 0.40]
Ypos_miranda = [0.35, 0.32, 0.30]

for i in range(len(Significance['targetscan'])):
    ax1.text(Xpos[i], Ypos_targetscan[i], Significance['targetscan'][i], horizontalalignment = 'center', verticalalignment = 'center', color = 'black', size = 8)
    ax2.text(Xpos[i], Ypos_miranda[i], Significance['miranda'][i], horizontalalignment = 'center', verticalalignment = 'center', color = 'black', size = 8)

# add legend relative to ax1 using ax1 coordinates
C = mpatches.Patch(facecolor = '#a6cee3', edgecolor = 'black', linewidth = 1, label= 'CNV')
N = mpatches.Patch(facecolor = '#b2df8a', edgecolor = 'black', linewidth = 1, label= 'non-CNV')
ax1.legend(handles = [C, N], loc = (-0.25, 1.2), fontsize = 8, frameon = False, ncol = 2)

# make sure subplots do not overlap
plt.tight_layout()

# build outputfile with arguments
outputfile = 'PlotTargetsGoldMap_' + CNVFilter.capitalize() + '_' + chromos + '_' + cnv_length
print(outputfile)

# save figure
fig.savefig(outputfile + '.eps', bbox_inches = 'tight')    