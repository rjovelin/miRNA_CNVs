# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 19:16:32 2015

@author: Richard
"""

import sys

# usage python3 merge_chicken_CNV_SuppData.py [long_CNVs/all_CNVs] outputfile


# merge all CNV genes into a text file

# get option to consider all CNVs or CNV greater than 1 KB (all or long)
long_CNV = sys.argv[1] 

# get CNV genes from Jia et al 2013 Animal Genetics Supp data
infile = open('Jia_AnGen_2013_age12009-sup-0001-TableS1.txt', 'r')
# create a set of genes
JiaSuppData = set()

# skip headers
infile.readline()
infile.readline()
# loop over file
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split()
        # get CNV length
        CNV_length = int(line[1].replace(',', ''))
        # get the CNV genes
        CNV_genes = line[-1].split(',')
        # check CNV length threshold
        if long_CNV == 'long_CNVs' and CNV_length >= 1000:
            # add CNV genes to set
            for gene in CNV_genes:
                JiaSuppData.add(gene)
        elif long_CNV == 'all_CNVs':
            for gene in CNV_genes:
                JiaSuppData.add(gene)
# close file
infile.close()
print('Jia Supp', len(JiaSuppData))


# get CNV_genes from Yi et al 2013 BMC Genomics 
infile = open('Yi_BMCGen_2014_1471-2164-15-962-s8.txt', 'r')
# create set of genes
YiSuppData = set()

# skip header
infile.readline()
# loop over file
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split()
        # get start, and position 0-bsed
        start = int(line[2]) - 1
        end = int(line[3])
        # get CNV length
        CNV_length = end - start
        # get CNV gene
        CNV_gene = line[6]
        # check CNV length threshold
        if long_CNV == 'long_CNVs' and CNV_length >= 1000:
            # add gene to set
            YiSuppData.add(CNV_gene)
        elif long_CNV == 'all_CNVs':
            # add gene to set
            YiSuppData.add(CNV_gene)
            
# close file
infile.close()
print('Yi Supp', len(YiSuppData))


# get CNV genes from han et al 2014 BMC Genomics           
infile = open('Han_BMCGen_2014_1471-2164-15-934-s3.txt', 'r')
# create set of genes
HanSuppData = set()

# skip headers
infile.readline()
infile.readline()
# loop over file
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split()
        # get CNV start and end position 0-based
        start = int(line[2]) - 1
        end = int(line[3])
        # get CNV length
        CNV_length = end -  start
        # check if gene is defined
        if len(line) > 5:
            # grab CNV gene name
            CNV_gene = line[6]
        # check CNV length threshold
        if long_CNV == 'long_CNVs' and CNV_length > 1000:
            # add gene to set
            HanSuppData.add(CNV_gene)
        elif long_CNV == 'all_CNVs':
            HanSuppData.add(CNV_gene)

# close file
infile.close()
print('Han Supp', len(HanSuppData))


# get CNV genes from Crooijmans et al 2013 BMC Genomics
infile = open('Crooijmans_BMCGen_2013_1471-2164-14-398-s5.txt', 'r')
# create set o genes
CrooijmansSuppData = set()
# skip headers
infile.readline()
infile.readline()
# loop over file
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split()
        # get CNV start and end positions 0-based
        start = int(line[3]) -1
        end = int(line[4])
        # get CNV length
        CNV_length = end - start
        # check if gene is defined
        if len(line) > 5:
            # get CNV gene
            CNV_gene = line[5]
            # remove _CHICK from gene name
            if '_CHICK' in CNV_gene:
                CNV_gene = CNV_gene[:CNV_gene.index('_CHICK')]
            # check CNV length threshold
            if long_CNV == 'long_CNVs' and CNV_length > 1000:
                # add gene to set
                CrooijmansSuppData.add(CNV_gene)
            elif long_CNV == 'all_CNVs':
                CrooijmansSuppData.add(CNV_gene)

# close file
infile.close()
print('Crooijmans Supp', len(CrooijmansSuppData))


# get CNV genes from Fan et al 2013 GBE Table S7
infile = open('Fan_GBE_2013_Supplementary_table_S7.txt', 'r')
# create a set of CNV genes
FanSuppData = set()
# skip header
infile.readline()
# loop over file
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split()
        # get CNV start and end position
        start = int(line[1]) -1
        end = int(line[2])
        # get CNV length
        CNV_length = end - start
        # check if gene is defined
        if len(line) == 10:
            # get CNV gene name
            CNV_gene = line[4]
            # trim name from _CHICK
            if '_CHICK' in CNV_gene:
                CNV_gene = CNV_gene[: CNV_gene.index('_CHICK')]
            # check CNV length threshold
            if long_CNV == 'long_CNVs' and CNV_length > 1000:
                # add gene to set
                FanSuppData.add(CNV_gene)
            elif long_CNV == 'all_CNVs':
                FanSuppData.add(CNV_gene)

# close file
infile.close()
print('Fan Supp', len(FanSuppData))


# create a set of CNV genes by merging all Supp data
merged_CNV_genes = set()

# create a list with set of CNV genes
SuppData = [JiaSuppData, YiSuppData, HanSuppData, CrooijmansSuppData, FanSuppData]

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