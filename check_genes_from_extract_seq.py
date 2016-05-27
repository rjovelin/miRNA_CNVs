# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 12:04:04 2015

@author: RJovelin
"""

# usage check_genes_from_extract_seq.py [True/False]


# verify that functions extracting UTR and CDS sequences get the same ganes


from CNV_miRNAs import *
import os
import sys

# get the option to keep genes on all chromos (False) or only on assembled 
# nuclear chromosomes only from the command

keep_valid_chromos = sys.argv[1]
if keep_valid_chromos == 'True':
    keep_valid_chromos = True
elif keep_valid_chromos == 'False':
    keep_valid_chromos = False
print(keep_valid_chromos)

# make a dictionary of species names : species code
species_names = {'H_sapiens': 'Hsa',  'P_troglodytes': 'Ptr', 'M_mulatta': 'Mmul',
                 'M_musculus': 'Mmus', 'R_norvegicus': 'Rno', 'B_taurus': 'Bta',
                 'C_familiaris': 'Cfa', 'G_gallus': 'Gga'}

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
    
    # extract the 3'UTR sequences from genome {gene name : UTR sequence}
    UTR_seq3 = extract_UTR_sequences(GFF, fasta, valid_chromos, '3UTR', keep_valid_chromos)
    print(species, '3UTR', len(UTR_seq3))    
    UTR_seq5 = extract_UTR_sequences(GFF, fasta, valid_chromos, '5UTR', keep_valid_chromos)
    print(species, '5UTR', len(UTR_seq5))
    CDS_seq = extract_CDS_sequences(GFF, fasta, valid_chromos, keep_valid_chromos)
    print(species, 'CDS', len(CDS_seq))
    # get gene names
    UTR_names3 = set(i for i in UTR_seq3)
    UTR_names5 = set(i for i in UTR_seq5)
    CDS_names = set(i for i in CDS_seq)
    print(species, '3UTR-5UTR', UTR_names3 == UTR_names5)
    print(species, '3UTR-CDS', UTR_names3 == CDS_names)
    print(species, '5UTR-CDS', UTR_names5 == CDS_names)






     