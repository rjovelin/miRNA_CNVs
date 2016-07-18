# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 00:04:04 2015

@author: Richard
"""

# usage: python3 merge_cow_CNV_SuppData.py [long_CNVs/all_CNVs] outputfile


import sys

# merge all CNV genes into a text file

# get option to consider all CNVs or CNV greater than 1 KB (all or long)
long_CNV = sys.argv[1] 


# get CNV genes from Brichart et al 2012 Genome Research Supp data
infile = open('Brickhart_GenRes_2012_Table_S7_CNVR_RefSeq.txt', 'r')
# create a set of genes
BrichartSuppData = set()

# skip headers
infile.readline()
# loop over file
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split()
        # get the CNV genes
        CNV_gene = line[0]
        BrichartSuppData.add(CNV_gene)
        
# close file
infile.close()
print('Brichart Supp', len(BrichartSuppData))


# get CNV_genes from Cicconardi et al 2014 BMC Genomics 
infile = open('Cicconardi_BMCGen_2014_1471-2164-14-124-s4.tsv', 'r')
# create set of genes
CicconardiSuppData = set()
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
        # get CNV gene
        CNV_gene = line[-1]
        # check that gene is defined
        if CNV_gene != '-':
            # remove BoVIN from gene name
            if '_BOVIN' in CNV_gene:
                CNV_gene = CNV_gene[:CNV_gene.index('_BOVIN')]
            # check CNV length threshold
            if long_CNV == 'long_CNVs' and CNV_length >= 1000:
                # add gene to set
                CicconardiSuppData.add(CNV_gene)
            elif long_CNV == 'all_CNVs':
                # add gene to set
                CicconardiSuppData.add(CNV_gene)
            
# close file
infile.close()
print('Cicconardi Supp', len(CicconardiSuppData))


# get CNV genes from Hou et al 2011 BMC Genomics           
infile = open('Hou_BMCGenom_2011_1471-2164-12-127-s2.txt', 'r')
# create set of genes
HouSuppData = set()

# skip headers
infile.readline()
infile.readline()
# loop over file
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split('\t')
        # get CNV length
        CNV_length = int(line[6])
        # get CNV genes
        if ';' in line[-2]:
            # multiple genes
            CNV_genes = line[-2].split(';')
            # extract the gene name
            for i in CNV_genes:
                i = i.split()
                gene = i[1]
                # check CNV length threshold
                if long_CNV == 'long_CNVs' and CNV_length >= 1000:
                    HouSuppData.add(gene)
                elif long_CNV == 'all_CNVs':
                    HouSuppData.add(gene)                
        else:
            CNV_genes = line[-2].split()
            gene = CNV_genes[1]
            # check CNV length threshold
            if long_CNV == 'long_CNVs' and CNV_length >= 1000:
                HouSuppData.add(gene)
            elif long_CNV == 'all_CNVs':
                HouSuppData.add(gene)  
            
# cloose file
infile.close()
print('Hou Supp', len(HouSuppData))


# create a set of CNV genes by merging all Supp data
merged_CNV_genes = set()

# create a list with set of CNV genes
SuppData = [BrichartSuppData, CicconardiSuppData, HouSuppData]

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