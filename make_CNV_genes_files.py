# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 16:37:54 2015

@author: RJovelin
"""




# usage make_CNV_genes_files.py [True/False] [long_CNVs/all_CNVs]


# generate files with CNV status for all genes used to predict miRNA target sites


from CNV_miRNAs import *
import os
import sys

# get the option to keep genes on all chromos (False) or only on assembled 
# nuclear chromosomes only (True) from the command
keep_valid_chromos = sys.argv[1]
if keep_valid_chromos == 'True':
    keep_valid_chromos = True
    chromos = 'valid_chromos'
elif keep_valid_chromos == 'False':
    keep_valid_chromos = False
    chromos = 'all_chromos'
print(keep_valid_chromos)

# get the option to call a CNV if CNV length > 1 Kb (long_CNVs)
# or to include all CNVs regardless of length (all_CNVs)
long_CNV = sys.argv[2]
print(long_CNV)
if long_CNV == 'all_CNVs':
    cnv_length = 'CNV_all_length'
elif long_CNV == 'long_CNVs':
    cnv_length = 'CNV_greater_1Kb'


# make a dictionary of species  : directory
species_dir = {'H_sapiens': './CNV_genes_published_data/Human/',
                 'P_troglodytes': './CNV_genes_published_data/Chimp/',
                 'M_mulatta': './CNV_genes_published_data/Macaque/',
                 'M_musculus': './CNV_genes_published_data/Mouse/',
                 'R_norvegicus': './CNV_genes_published_data/Rat/',
                 'B_taurus': './CNV_genes_published_data/Cow/',
                 'C_familiaris': './CNV_genes_published_data/Dog/',
                 'G_gallus': './CNV_genes_published_data/Chicken/'}

# make a dictionary of species : common name
species_names = {'H_sapiens': 'human', 'P_troglodytes': 'chimp',
                 'M_mulatta': 'macaque', 'M_musculus': 'mouse',
                 'R_norvegicus': 'rat', 'B_taurus': 'cow',
                 'C_familiaris': 'dog', 'G_gallus': 'chicken'}

# loop over species name
for species in species_names:
    print(species)
    # get GFF file
    GFF = species + '.gff3'
    print(GFF)
    # get genome file
    fasta = species + '_genome.txt'
    print(fasta)
    # get the file with valid chromos
    valid_chromos = species + '_valid_chromos.txt'
    print(valid_chromos)
    
    # check if species is human
    if species == 'H_sapiens':
        # generate CNV status files for each version of the DGV
        # make a list of files with cnv genes
        cnv_files = ['human_CNV_all_length_GRCh37_2014.txt',
                     'human_CNV_all_length_GRCh37_2013-05.txt',
                     'human_CNV_all_length_GRCh37_2013-07.txt',
                     'human_CNV_all_length_GRCh37_2015.txt']
                     
        # generate CNV status files for each version of the DGV
        for filename in cnv_files:
            # get outputfile
            if '2013' in filename:
                outputfile = species + '_GRCh37_' + filename[filename.index('2013'): filename.index('.txt')] + '_' + cnv_length + '_' + chromos + '.txt'
            else:
                outputfile = species + '_GRCh37_' + filename[filename.rindex('_') + 1: filename.index('.txt')] + '_' + cnv_length + '_' + chromos + '.txt'
            print(outputfile)
            
            # get CNV status
            CNV_status = get_genes_CNV_status(GFF, fasta, valid_chromos, keep_valid_chromos, species_dir[species] + filename)
            print(filename, len(CNV_status))
            
            # open file for writing
            newfile = open(outputfile, 'w')
            # write header
            newfile.write('gene\tCNV_status\n')
            # loop over gene in CNV_status
            for gene in CNV_status:
                # write gene and status to file
                newfile.write(gene + '\t' + CNV_status[gene] + '\n')
            newfile.close()
            
            
    else:
        # get the CNV file
        CNV_file = species_names[species] + '_' + cnv_length + '.txt'
        print(species, CNV_file)
        
        # get CNV status for all genes
        CNV_status = get_genes_CNV_status(GFF, fasta, valid_chromos, keep_valid_chromos, species_dir[species] + CNV_file)
        print(species, len(CNV_status))
        
        # get outputfile
        outputfile = species + '_' + cnv_length + '_' + chromos + '.txt'
        print(species, outputfile)
        
        # open file for writing
        newfile = open(outputfile, 'w')
        # write header
        newfile.write('gene\tCNV_status\n')
        
        # loop over gene in CNV_status
        for gene in CNV_status:
            # write gene and status to file
            newfile.write(gene + '\t' + CNV_status[gene] + '\n')
            
        # close file after writing
        newfile.close()
        print('done sorting {0} genes according to CNV status'.format(species))        

