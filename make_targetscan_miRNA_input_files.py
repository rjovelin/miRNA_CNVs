# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 15:59:27 2015

@author: RJovelin
"""

from CNV_miRNAs import *
import os

# make a dictionary of species names : species code
species_names = {'H_sapiens': 'Hsa',  'P_troglodytes': 'Ptr', 'M_mulatta': 'Mmul',
                 'M_musculus': 'Mmus', 'R_norvegicus': 'Rno', 'B_taurus': 'Bta',
                 'C_familiaris': 'Cfa', 'G_gallus': 'Gga'}

# loop over species name
for species in species_names:
    print(species, species_names[species])
    # get mature file
    mature = species + '_mature.txt'
    print(mature)
    
    outputfile = species + '_miRFam_targetscan.txt'    
    print(outputfile)
    
    # make mIRNA family Targetscan input file
    make_targetscan_mirfam_input_file(mature, species_names[species], outputfile)
    print('done making {0} Targetscan miRNA input file'.format(outputfile))    
    
    