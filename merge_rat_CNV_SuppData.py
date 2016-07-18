# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 11:16:34 2015

@author: RJovelin
"""

# usage: python3 merge_rat_CNV_SuppData.py [long_CNVs/all_CNVs] outputfile


import sys

# merge all CNV genes into a text file

# get option to consider all CNVs or CNV greater than 1 KB (all or long)
long_CNV = sys.argv[1] 


# make a dictionary with Rat gene coordinates {chromo: gene: [start:end]}
Rat_genes = {}
# open annotation file for reading
infile = open('R_norvegicus_rn4_RefSeq_Genes.txt', 'r')
# skip header
infile.readline()
# loop ober file
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split('\t')
        # get chromo (eg. chr1)
        chromo = line[2]
        # get gene name
        gene = line[12]
        # get gene start and end position (already in 0-based in file)
        start = int(line[4])
        end = int(line[5])
        # check if chromo in outter dict
        if chromo not in Rat_genes:
            # initialize inner dict
            Rat_genes[chromo] = {}
        # populate dict with gene : coordinate pair
        Rat_genes[chromo][gene] = [start, end]
        
# close file
infile.close()
print(len(Rat_genes))

# create a dict with CNV coordinates from Guryev et al 2008 Nature Genetics Table S3
# {chromo: set(position of CNV on chromo)}
CNVS3_coord = {}

# open file for reading
infile = open('Guryev_NatGenet_2008_tableS3.txt', 'r')
# skip header
infile.readline()
# loop over file
for line in infile:
    # remove left and right white space
    line = line.strip()
    if line != '':
        line = line.split()
        # get chromo
        chromo = 'chr' + line[0]
        # get CNV start and end position 0-based
        # remove comma in positions
        start = int(line[1].replace(',', '')) -1
        end = int(line[2].replace(',', ''))
        # check if chromo in dict
        if chromo not in CNVS3_coord:
            # add key, initialize set
            CNVS3_coord[chromo] = set()
        # check CNV length threshold
        if long_CNV == 'long_CNVs' and end - start > 1000:
            # add CNV positions to chromo
            for i in range(start, end):
                CNVS3_coord[chromo].add(i)
        elif long_CNV == 'all_CNVs':
            # add CNV positions to chromo
            for i in range(start, end):
                CNVS3_coord[chromo].add(i)
            
# close file
infile.close()

# remove chromo without positions
to_remove = []
for chromo in CNVS3_coord:
    if len(CNVS3_coord[chromo]) == 0:
        to_remove.append(chromo)
if len(to_remove) != 0:
    for chromo in to_remove:
        del CNVS3_coord[chromo]        

print(len(CNVS3_coord))
# print chromo and number of sites
for chromo in CNVS3_coord:
    print('CNVS3_coord', chromo, len(CNVS3_coord[chromo]))
    

# create a dict with CNV coordinates from Guryev et al 2008 Nature Genetics Table S4
# {chromo: set(position of CNVs on chromo)}
CNVS4_coord = {}

# open file for reading
infile = open('Guryev_NatGenet_2008_tableS4.txt', 'r')
# skip header
infile.readline()
# loop over file
for line in infile:
    # remove left and right white space
    line = line.strip()
    if line != '':
        line = line.split()
        # get chromo
        chromo = 'chr' + line[0]
        # get CNV start and end position 0-based
        # remove comma in positions
        start = int(line[1].replace(',', '')) -1
        end = int(line[2].replace(',', ''))
        # check if chromo in dict
        if chromo not in CNVS4_coord:
            # add key, initialize set
            CNVS4_coord[chromo] = set()
        # check CNV length threshold
        if long_CNV == 'long_CNVs' and end - start > 1000:
            # add CNV positions to chromo
            for i in range(start, end):
                CNVS4_coord[chromo].add(i)
        elif long_CNV == 'all_CNVs':
            # add CNV positions to chromo
            for i in range(start, end):
                CNVS4_coord[chromo].add(i)
        
# close file
infile.close()

# remove chromo without positions
to_remove = []
for chromo in CNVS4_coord:
    if len(CNVS4_coord[chromo]) == 0:
        to_remove.append(chromo)
if len(to_remove) != 0:
    for chromo in to_remove:
        del CNVS4_coord[chromo]        

print(len(CNVS4_coord))
# print chromo and number of sites
for chromo in CNVS4_coord:
    print('CNVS4_coord', chromo, len(CNVS4_coord[chromo]))

# create sets of CNV genes 
CNVgeneS3, CNVgeneS4 = set(), set()


# loop over chromo in Rat genes
for chromo in Rat_genes:
    print(chromo)
    # chek that chromo in CNV coord
    if chromo in CNVS3_coord:
        # loop over genes on that chromo
        for gene in Rat_genes[chromo]:
            # get gene coordinates
            gene_pos = set(range(Rat_genes[chromo][gene][0], Rat_genes[chromo][gene][1]))
            # compare gene position with CNV positions
            if len(gene_pos.intersection(CNVS3_coord[chromo])) != 0:
                CNVgeneS3.add(gene)

print('CNV genes from Table S3', len(CNVgeneS3))


# loop over chromo in Rat_genes
for chromo in Rat_genes:
    print(chromo)
    # check that chromo in CNV_coord:
    if chromo in CNVS4_coord:
        # loop over genes on that chromo
        for gene in Rat_genes[chromo]:
            # get gene coordinates
            gene_pos = set(range(Rat_genes[chromo][gene][0], Rat_genes[chromo][gene][1]))
            # compare gene position with CNV positions
            if len(gene_pos.intersection(CNVS4_coord[chromo])) != 0:
                CNVgeneS4.add(gene)

print('CNV genes from Table S4', len(CNVgeneS4))

# create a set of genes by merging the 2 sets of CNV genes
merged_CNV_genes = CNVgeneS3.union(CNVgeneS4)
# make a list from set
merged_CNV_genes = list(merged_CNV_genes)
# sort list
merged_CNV_genes.sort()
print(len(merged_CNV_genes))

# open file for writing
# get outputfile from the command
outputfile = sys.argv[2]

newfile = open(outputfile, 'w')
# loop over list of genes
for gene in merged_CNV_genes:
    # write content to file
    newfile.write(gene + '\n')

# close file
newfile.close()