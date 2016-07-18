# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 14:40:43 2015

@author: RJovelin
"""

# usage: python3 merge_macaque_CNV_SuppData.py [long_CNVs/all_CNVs] outputfile


import sys

# merge all CNV genes into a text file

# get option to consider all CNVs or CNV greater than 1 KB (all or long)
long_CNV = sys.argv[1] 

# make a dictionary with Macaque gene coordinates {chromo: gene: [start:end]}
rhesus_genes = {}
# open annotation file for reading
infile = open('M_mulatta_rheMac2_RefSeq_Annot.txt', 'r')
# skip header
infile.readline()
# loop ober file
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split('\t')
        assert len(line) == 16, 'line length is not constant, inaccurate gene position'
        # get chromo (eg. chr1)
        chromo = line[2]
        # get gene name
        gene = line[12]
        # get gene start and end position (already in 0-based in file)
        start = int(line[4])
        end = int(line[5])
        # check if chromo in outter dict
        if chromo not in rhesus_genes:
            # initialize inner dict
            rhesus_genes[chromo] = {}
        # populate dict with gene : coordinate pair
        rhesus_genes[chromo][gene] = [start, end]
        
# close file
infile.close()
print('rhesus chromo', len(rhesus_genes))


# create a dict with CNV coordinates from Gokcumen et al 2013 PNAS CNVs
# {chromo: set(positions of CNV on that chromo)
rhesus_CNV_coord = {}

# open file with deletions for reading
infile = open('Gokcumen_Macaque_Del.txt', 'r')
# skip header
infile.readline()
# loop over file
for line in infile:
    # remove right white space
    line = line.rstrip()
    if line != '':
        line = line.split('\t')
        # get chromo
        chromo = line[2]
        # get CNV start and end positions
        start = int(line[3]) - 1
        end = int(line[4])
        # check if chromo in dict
        if chromo not in rhesus_CNV_coord:
            # initialise set value
            rhesus_CNV_coord[chromo] = set()
        # check CNV length threshold
        if long_CNV == 'long_CNVs' and end - start > 1000:
            # add CNV positions to chromo
            for i in range(start, end):
                rhesus_CNV_coord[chromo].add(i)
        elif long_CNV == 'all_CNVs':
            # add CNV_position to chromo
            for i in range(start, end):
                rhesus_CNV_coord[chromo].add(i)
                
# close file
infile.close()
print('rhesus del CNV', len(rhesus_CNV_coord))


# open file with duplications for reading
infile = open('Gokcumen_Macaque_Dupli.txt', 'r')
# skip header
infile.readline()
# loop over file
for line in infile:
    # remove left and right white space
    line = line.rstrip()
    if line != '':
        line = line.split('\t')
        # get chromo
        chromo = line[2]
        # get CNV start and end positions
        start = int(line[3]) - 1
        end = int(line[4])
        # check if chromo in dict
        if chromo not in rhesus_CNV_coord:
            # initialise set value
            rhesus_CNV_coord[chromo] = set()
        # check CNV length threshold
        if long_CNV == 'long_CNVs' and end - start > 1000:
            # add CNV positions to chromo
            for i in range(start, end):
                rhesus_CNV_coord[chromo].add(i)
        elif long_CNV == 'all_CNVs':
            # add CNV_position to chromo
            for i in range(start, end):
                rhesus_CNV_coord[chromo].add(i)        
        
# close file
infile.close()
print('rhesus dupli and del CNV', len(rhesus_CNV_coord))


# add CNV coordinates for new MEIs

# open file with new MEIs for reading
infile = open('Gokcumen_Macaque_newMEI.txt', 'r')
# skip header
infile.readline()
# loop over file
for line in infile:
    # remove left and right white space
    line = line.rstrip()
    if line != '':
        line = line.split('\t')
        # get chromo
        chromo = line[0]
        # get CNV start and end positions
        start = int(line[1]) -1
        end = int(line[2])
        # check if chromo in dict
        if chromo not in rhesus_CNV_coord:
            # initialise set value
            rhesus_CNV_coord[chromo] = set()
        # check CNV length threshold
        if long_CNV == 'long_CNVs' and end - start > 1000:
            # add CNV positions to chromo
            for i in range(start, end):
                rhesus_CNV_coord[chromo].add(i)
        elif long_CNV == 'all_CNVs':
            # add CNV_position to chromo
            for i in range(start, end):
                rhesus_CNV_coord[chromo].add(i)        
        
# close file
infile.close()
print('rhesus all CNV', len(rhesus_CNV_coord))


# create a set of CNV genes
rhesus_CNV_genes = set()

# loop over chromo in rhesus genes
for chromo in rhesus_genes:
    print(chromo)
    # check that chromo in CNV coord
    if chromo in rhesus_CNV_coord:
        # loop over genes on chromo
        for gene in rhesus_genes[chromo]:
            # get the positions
            gene_pos = set(range(rhesus_genes[chromo][gene][0], rhesus_genes[chromo][gene][1]))
            # compare gene position with CNV positions
            if len(gene_pos.intersection(rhesus_CNV_coord[chromo])) != 0:
                rhesus_CNV_genes.add(gene)
                
                
print('N rhesus CNV genes', len(rhesus_CNV_genes))

# make a list from set
rhesus_CNV_genes = list(rhesus_CNV_genes)
# sort list
rhesus_CNV_genes.sort()

# open file for writing
# get outputfile from the command
outputfile = sys.argv[2]

newfile = open(outputfile, 'w')
# loop over list of genes
for gene in rhesus_CNV_genes:
    # write content to file
    newfile.write(gene + '\n')

# close file
newfile.close()