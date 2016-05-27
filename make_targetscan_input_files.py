# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 14:21:35 2015

@author: RJovelin
"""


# usage python3 make_targetscan_input_files.py [True/False]


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


# make a dictionary of species names : species code
species_names = {'H_sapiens': 'Hsa',  'P_troglodytes': 'Ptr', 'M_mulatta': 'Mmul',
                 'M_musculus': 'Mmus', 'R_norvegicus': 'Rno', 'B_taurus': 'Bta',
                 'C_familiaris': 'Cfa', 'G_gallus': 'Gga'}

# make a list of gene regions
gene_regions = ['5UTR', '3UTR', 'CDS']

# loop over gene regions, make input file for each gene region
for region in gene_regions:
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
            
        # check if keep genes on all chromos or on assembled nuclear chromos only
        if keep_valid_chromos == True:
            outputfile = species + '_' + region + '_valid_chromos_targetscan.txt'
        elif keep_valid_chromos == False:
            outputfile = species + '_' + region + '_all_chromos_targetscan.txt'
        print(outputfile)        
        
        # generate Targetscan sequence input file
        make_targetscan_seq_input_file(GFF, fasta, valid_chromos, keep_valid_chromos, region, species_names[species], outputfile)
        print('done making {0} TargetScan input file'.format(outputfile))        
        
