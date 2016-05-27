# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 16:31:18 2015

@author: RJovelin
"""

# usage python3 make_genes_3UTR_length_files.py [True/False]


# sort genes into short 3' UTR or long 3' UTR


from CNV_miRNAs import *
import os
import sys

# get the option to keep genes on all chromos (False) or only on assembled 
# nuclear chromosomes only from the command

keep_valid_chromos = sys.argv[1]
if keep_valid_chromos == 'True':
    keep_valid_chromos = True
    chromos = 'valid_chromos'
elif keep_valid_chromos == 'False':
    keep_valid_chromos = False
    chromos = 'all_chromos'
print(keep_valid_chromos)

# make a dictionary of species names : species code
species_names = {'H_sapiens': 'Hsa',  'P_troglodytes': 'Ptr', 'M_mulatta': 'Mmul',
                 'M_musculus': 'Mmus', 'R_norvegicus': 'Rno', 'B_taurus': 'Bta',
                 'C_familiaris': 'Cfa', 'G_gallus': 'Gga'}

# loop over species name
for species in species_names:
    print(species, species_names[species])
    # get GFF file
    GFF = species + '.gff3'
    print(GFF)
    # get genome file
    fasta = species + '_genome.txt'
    print(fasta)
       
    # get the file with valid chromos
    valid_chromos = species + '_valid_chromos.txt'
    print(valid_chromos)    
    
    # get outputfile
    outputfile = species + '_3UTR_length_' + chromos + '.txt'
    print(outputfile)        
        
    # extract the 3'UTR sequences from genome {gene name : UTR sequence}
    UTR_seq = extract_UTR_sequences(GFF, fasta, valid_chromos, '3UTR', keep_valid_chromos)
    print(species, len(UTR_seq))    
    
    # open file for writing
    newfile = open(outputfile, 'w')
    # write header
    newfile.write('gene\t3UTR_type\t3UTR_length\n')

    # loop over gene in UTR_seq:
    for gene in UTR_seq:
        # check if UTR length < 7bp (minimum requirement for a miRNA target site)
        if len(UTR_seq[gene]) < 7:
            # gene has short 3' UTR
            newfile.write('\t'.join([gene, 'short', str(len(UTR_seq[gene]))]) + '\n')
        else:
            # gene is short
            newfile.write('\t'.join([gene, 'long', str(len(UTR_seq[gene]))]) + '\n')
            
    # close file after writing
    newfile.close()
    print('done sorting {0} genes according to UTR length'.format(species))        

