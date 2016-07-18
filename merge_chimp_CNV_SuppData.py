# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 19:50:10 2015

@author: Richard
"""

# usage: python3 merge_chimp_CNV_SuppData [long_CNVs/all_CNVs] outputfile


import sys

# merge all CNV genes into a text file

# get option to consider all CNVs or CNV greater than 1 KB (all or long)
long_CNV = sys.argv[1] 


# make a dictionary with chimp ensembl genes coordinates {chromo: gene: [start:end]}
# annotations from panTro2 have been Lifted-Over (using USCS Lift-Over tool) to panTro3 coordinates

chimp_genes = {}
# open annotation file for reading
infile = open('hglft_pantro2_to_pantro3.txt', 'r')
# skip header
infile.readline()
# loop ober file
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split('\t')
        # get chromo
        chromo = line[0]
        # get gene name
        gene = line[3]
        # get gene start and end position (already in 0-based in file)
        start = int(line[1])
        end = int(line[2])
        # check if chromo in outter dict
        if chromo not in chimp_genes:
            # initialize inner dict
            chimp_genes[chromo] = {}
        # populate dict with gene : coordinate pair
        chimp_genes[chromo][gene] = [start, end]
        
# close file
infile.close()
print('chimp chromo', len(chimp_genes))

# count the number of genes
total = 0
for chromo in chimp_genes:
    for gene in chimp_genes[chromo]:
        total += 1
print('N chimp genes', total)


# create a dict with CNV coordinates from Gokcumen et al 2013 PNAS CNVs
# {chromo: set(position of CNVs on that chromo)
chimp_CNV_coord = {}

# open file with deletions for reading
infile = open('Gokcumen_Chimp_Del.txt', 'r')
# skip header
infile.readline()
# loop over file
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split('\t')
        # get chromo
        chromo = line[2]
        # get CNV start and end positions
        start = int(line[3]) - 1
        end = int(line[4])
        # check if chromo in dict
        if chromo not in chimp_CNV_coord:
            # initialise set value
            chimp_CNV_coord[chromo] = set()
        # check CNV length threshold
        if long_CNV == 'long_CNVs' and end - start > 1000:
            # add CNV positions to chromo
            for i in range(start, end):
                chimp_CNV_coord[chromo].add(i)
        elif long_CNV == 'all_CNVs':
            # add CNV positions to chromo
            for i in range(start, end):
                chimp_CNV_coord[chromo].add(i)

# close file
infile.close()
print('chimp del CNV', len(chimp_CNV_coord))



# open file with duplications for reading
infile = open('Gokcumen_Chimp_Dupli.txt', 'r')
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
        if chromo not in chimp_CNV_coord:
            # initialize set value
            chimp_CNV_coord[chromo] = set()
        # check CNV length threshold
        if long_CNV == 'long_CNVs' and end - start > 1000:
            # add CNV positions to chromo
            for i in range(start, end):
                chimp_CNV_coord[chromo].add(i)
        elif long_CNV == 'all_CNVs':
            # add CNV positions to chromo
            for i in range(start, end):
                chimp_CNV_coord[chromo].add(i)

# close file
infile.close()
print('chimp dupli and del CNV', len(chimp_CNV_coord))


# add CNV coordinates for new MEIs

# open file with new MEIs for reading
infile = open('Gokcumen_Chimp_newMEI.txt', 'r')
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
        if chromo not in chimp_CNV_coord:
            chimp_CNV_coord[chromo] = set()
        # check CNV length threshold
        if long_CNV == 'long_CNVs' and end - start > 1000:
            # add CNV positions to chromo
            for i in range(start, end):
                chimp_CNV_coord[chromo].add(i)
        elif long_CNV == 'all_CNVs':
            # add CNV positions to chromo
            chimp_CNV_coord[chromo].add(i)

# close file
infile.close()
print('chimp all CNV', len(chimp_CNV_coord))

# remove chromo without positions
to_remove = []
for chromo in chimp_CNV_coord:
    if len(chimp_CNV_coord[chromo]) == 0:
        to_remove.append(chromo)
if len(to_remove) != 0:
    for chromo in to_remove:
        del chimp_CNV_coord[chromo] 
print('chimp all CNV', len(chimp_CNV_coord))


# create a set of CNV genes
chimp_CNV_genes = set()


# loop over chromo in chimp_genes
for chromo in chimp_genes:
    print(chromo)
    # check that chromo in CNV coord
    if chromo in chimp_CNV_coord:
        # loop over genes on that chromo
        for gene in chimp_genes[chromo]:
            # get gene positions
            gene_pos = set(range(chimp_genes[chromo][gene][0], chimp_genes[chromo][gene][1]))
            # compare gene position with CNV positions
            if len(gene_pos.intersection(chimp_CNV_coord[chromo])) != 0:
                chimp_CNV_genes.add(gene)

                 
print('N chimp CNV genes', len(chimp_CNV_genes))


# make a dict ENSEMBL transcript : Refseq gene name
ensToName = {}
# open file for reading
infile = open('P_troglodytes_pantro2_ensToGeneName.txt', 'r')
# skip header
infile.readline()
# loop over file
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split()
        transcript = line[0]
        gene = line[1]
        ensToName[transcript] = gene
        
# close file after reading
infile.close()


# make a dict withENSEMNL gene name : [list of ENSEMBL transcripts]
ensGeneTotranscript = {}

# open file for reading
infile = open('P_troglodytes_pantro2_ensGtp.txt', 'r')
# skip header
infile.readline()
# loop over file
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split()
        # get gene
        gene = line[0]
        # get transcript
        transcript = line[1]
        # check if gene is already recorded
        if gene in ensGeneTotranscript:
            # add transcript to list
            ensGeneTotranscript[gene].append(transcript)
        else:
            # add gene : [transcript] pair
            ensGeneTotranscript[gene] = [transcript]
            
# close file
infile.close()

# create a set of CNV RefSeq gene names
CNV_gene_names = set()

# loop over CNV genes
for gene in chimp_CNV_genes:
    # get all transcript names corresponding to that gene
    for transcript in ensGeneTotranscript[gene]:
        # check if transcript has a RefSeq gene name
        if transcript in ensToName:
            # add gane name to CNV set
            CNV_gene_names.add(ensToName[transcript])
            

print('CNV gene name', len(CNV_gene_names))

# make a list from set
CNV_gene_names = list(CNV_gene_names)
# sort list
CNV_gene_names.sort()

# open file for writing
# get outputfile from the command
outputfile = sys.argv[2]

newfile = open(outputfile, 'w')
# loop over list of genes
for gene in CNV_gene_names:
    # write content to file
    newfile.write(gene + '\n')

# close file
newfile.close()






























