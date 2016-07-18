# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 10:22:26 2015

@author: RJovelin
"""

# usage: python3 merge_mouse_CNV_SuppData.py [long_CNVs/all_CNVs] outputfile

import sys

# merge all CNV genes into a text file

# get option to consider all CNVs or CNV greater than 1 KB (all or long)
long_CNV = sys.argv[1] 


# get CNV genes from She et al 2008 Nature Genetics Supp data
infile = open('She_TableS6.txt', 'r')
# create a set of genes
SheSuppData = set()

# skip headers
infile.readline()
# loop over file
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split('\t')
        # get CNV start and end positions 0-based
        start = int(line[2]) -1
        end = int(line[3])
        # get CNV length
        CNV_length = end - start
        # get the CNV genes
        CNV_genes = line[19].split(',')
        # get all genes
        for gene in CNV_genes:
            # trim white space
            gene = gene.strip()
            # check CNV length threshold
            if long_CNV == 'long_CNVs' and CNV_length >= 1000:
                # add gene to set
                SheSuppData.add(gene)
            elif long_CNV == 'all_CNVs':
                # add gene to set
                SheSuppData.add(gene)            
            
# close file
infile.close()
print('She Supp', len(SheSuppData))


# get CNV_genes from Locke et al 2015 BMC Genomics 
infile = open('Locke_BMCGen_2015_s12864-015-1713-z-s3.txt', 'r')
# create set of genes
LockeSuppData = set()
# skip header
infile.readline()
# loop over file
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split()
        # get start, and position 0-bsed
        start = int(line[1]) - 1
        end = int(line[2])
        # get CNV length
        CNV_length = end - start
        # get CNV genes
        CNV_genes = line[-1].split(',')
        # get all genes
        for gene in CNV_genes:
            # trim white space
            gene = gene.strip()
            # check CNV length threshold
            if long_CNV == 'long_CNVs' and CNV_length >= 1000:
                # add gene to set
                LockeSuppData.add(gene)
            elif long_CNV == 'all_CNVs':
                # add gene to set
                LockeSuppData.add(gene)
        
# close file
infile.close()
print('Locke Supp', len(LockeSuppData))


# get CNV genes from Pezer et al 2015 Genome Research           
infile = open('Pezer_GR_2015_Supp_Table_4.txt', 'r')
# create set of genes
PezerSuppData = set()

# skip headers
infile.readline()
infile.readline()
infile.readline()
# loop over file
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split('\t')
        # get coordinates
        position = line[0]
        # get CNV start and end position 0-based
        position = position[position.index(':') + 1:]
        position = position.split('-')
        start = int(position[0]) -1
        end = int(position[1])
        # get CNV gene
        CNV_gene = line[1]
        # check CNV length threshold
        if long_CNV == 'long_CNVs' and CNV_length >= 1000:
            # add gene to set
            PezerSuppData.add(CNV_gene)
        elif long_CNV == 'all_CNVs':
            # add gene to set
            PezerSuppData.add(CNV_gene)
            
# cloose file
infile.close()
print('Pezer Supp', len(PezerSuppData))


# create a set of CNV genes by merging all Supp data
merged_CNV_genes = set()

# create a list with set of CNV genes
SuppData = [SheSuppData, LockeSuppData, PezerSuppData]

# loop over list
for i in range(len(SuppData)):
    # loop over genes in set
    for gene in SuppData[i]:
        # add gene to merged CNV set
        merged_CNV_genes.add(gene)

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