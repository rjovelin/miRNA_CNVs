# -*- coding: utf-8 -*-
"""
Created on Wed Dec 23 14:44:01 2015

@author: RJovelin
"""

# create files to compare the number of mirna binding sites among species


from CNV_miRNAs import *
import os
import sys
import matplotlib.pyplot as plt

# usage python3 make_species_comp_files.py [options]
# [5UTR/CDS] use 5'UTR or CDS sites
# [True/False] use valid chromosomes (True) or all chromsomes (False)
# [long_CNVs/all_CNVs] use all CNVs or CNVs > 1Kb
# [CNV/not_CNV] compare CNV genes or non-CNV genes
# human_CNV_file choose the human CNV file

# get the region to consider to predict target sites [5UTr or CDS]
domain = sys.argv[1]
print(domain)

# get the option to keep genes on all chromos (False) or only on assembled 
# nuclear chromosomes only from the command
keep_valid_chromos = sys.argv[2]

# check if all chromos (including unplaced, unlocated, and MT) are used
# or if only valid chromos are used 
# note: only CNVs on valid chromos are reported in DGV, so if all chromos are
# used it may introduce a bias by calling non CNV genes genes that cannot be detected

if keep_valid_chromos == 'True':
    keep_valid_chromos = True
    chromos = 'valid_chromos'
elif keep_valid_chromos == 'False':
    keep_valid_chromos = False
    chromos = 'all_chromos'
print(keep_valid_chromos, chromos)


# get the type of CNVs to consider from the command
long_CNV = sys.argv[3]

if long_CNV == 'all_CNVs':
    cnv_length = 'CNV_all_length'
elif long_CNV == 'long_CNVs':
    cnv_length = 'CNV_greater_1Kb'
print(long_CNV, cnv_length)

# compare CNV genes or non-CNV genes
gene_type = sys.argv[4]

# get the human CNV file from the command, so that comparisons can be made 
# between all species and different version of the human DVG
human_CNV_file = sys.argv[5]

# get realease version of the DGV
release_version = human_CNV_file[human_CNV_file.index('GRCh37') : human_CNV_file.index('_CNV')]

# make a list of species names
species = ['H_sapiens', 'P_troglodytes', 'M_mulatta', 'M_musculus', 'B_taurus', 'G_gallus']

species_codes = {'H_sapiens': 'Hsa', 'P_troglodytes': 'Ptr', 'M_mulatta': 'Mml',
                 'M_musculus': 'Mmu', 'B_taurus': 'Bta', 'G_gallus':'Gga'}

# create a dict for human-species comparison with a list of list of sites for each species
# {Hsa_species_name : [[sites_Hsa], [sites_sp2]]}
species_data = {}

# loop over non-human species
for i in range(1, len(species)):
    # compare human to each other species
    print(species[0], species[i])
    # get genome files
    genome_sp1 = species[0] + '_genome.txt'
    genome_sp2 = species[i] + '_genome.txt'
    # get chromos files
    valid_chromos_sp1 = species[0] + '_valid_chromos.txt'
    valid_chromos_sp2 = species[i] + '_valid_chromos.txt'
    # get GFF annotation files
    GFF_sp1 = species[0] + '.gff3'
    GFF_sp2 = species[i] + '.gff3'
    print(genome_sp1, genome_sp2)
    print(valid_chromos_sp1, valid_chromos_sp2)
    print(GFF_sp1, GFF_sp2)

    # get targetscan sequence input file
    seq_input_sp1 = species[0] + '_' + domain + '_' + chromos + '_targetscan.txt'
    seq_input_sp2 = species[i] + '_' + domain + '_' + chromos + '_targetscan.txt'
    print(seq_input_sp1, seq_input_sp2)
        
    # get predictor output
    predicted_targets_sp1 = species[0] + '_' + domain + '_' + chromos + '_predicted_sites_targetscan.txt'
    predicted_targets_sp2 = species[i] + '_' + domain + '_' + chromos + '_predicted_sites_targetscan.txt'
    print(predicted_targets_sp1, predicted_targets_sp2)        
                
    # get UTR files
    UTR_sp1 = species[0] + '_3UTR_length_' + chromos + '.txt'      
    UTR_sp2 = species[i] + '_3UTR_length_' + chromos + '.txt'
    print(UTR_sp1, UTR_sp2)
         
    # get CNV files
    CNV_file_sp1 = human_CNV_file
    CNV_file_sp2 = species[i] + '_' + cnv_length + '_' + chromos + '.txt'
            
    # get mature mirna files
    mature_sp1 = species[0] + '_mature.txt'
    mature_sp2 = species[i] + '_mature.txt'        
    print(mature_sp1, mature_sp2)
        
    # get a set of shared seeds between species pair
    shared_seeds = find_conserved_mirna_families(mature_sp1, mature_sp2)
        
    # get CNV gene status
    CNV_status1 = sort_genes_CNV_status(CNV_file_sp1)
    CNV_status2 = sort_genes_CNV_status(CNV_file_sp2)
    print('CNV status', len(CNV_status1), len(CNV_status2))
        
    # parse predictor outputfile {gene: [N_sites, seq_length, N_sites/seq_length]}
    # consider only shared mirna families between species pairs
    # filter targets for shared mirnas
    targets_sp1 = parse_targetscan_output(seq_input_sp1, predicted_targets_sp1, 'conserved', shared_seeds, mature_sp1)
    targets_sp2 = parse_targetscan_output(seq_input_sp2, predicted_targets_sp2, 'conserved', shared_seeds, mature_sp2)
    print('targets', len(targets_sp1), len(targets_sp2))
        
    # sort genes by UTR length
    UTR_length_sp1 = sort_genes_3UTR_length(UTR_sp1)
    UTR_length_sp2 = sort_genes_3UTR_length(UTR_sp2)
    print('UTR length', len(UTR_length_sp1), len(UTR_length_sp2))
        
    # create list of normalized targets
    Sp1_sites, Sp2_sites =  [], []
        
    # loop over genes in targets sp1
    for gene in targets_sp1:
        # record only genes with short 3'UTR
        if UTR_length_sp1[gene] == 'short':
            # check if CNV or non-CNV genes should be recorded
            if CNV_status1[gene] == gene_type:
                # record normalized sites
                Sp1_sites.append(targets_sp1[gene][2])
                
    # loop over genes in targets sp2
    for gene in targets_sp2:
        # record only genes with short 3'UTR
        if UTR_length_sp2[gene] == 'short':
            # check if CNV or non-CNV genes should be recorded
            if CNV_status2[gene] == gene_type:
                Sp2_sites.append(targets_sp2[gene][2])
    
    print('targets', len(Sp1_sites), len(Sp2_sites))
    
    # populate dict
    comp_names = species_codes[species[0]] + '_' + species_codes[species[i]]
    print(comp_names)
    # copy lists to avoid modifying inner lists
    species_data[comp_names] = [Sp1_sites[:], Sp2_sites[:]]
    

for comp_names in species_data:
    print(comp_names, len(species_data[comp_names][0]), len(species_data[comp_names][1]))

# create a file to store the data
outputfile = 'Data_for_species_comparison_fig_' + domain + '_' + chromos + '_' + long_CNV + '_' + release_version + '_' + gene_type + '.txt' 
print(outputfile)
newfile = open(outputfile, 'w')
# write header
newfile.write('\t'.join(['Species1_Species2', 'Species', 'Sites']) + '\n')

# loop over dictionnary
for comp_names in species_data:
    # get species1 and species2
    species1, species2 = comp_names.split('_')
    # write comp_names, species_1 and corresponding sites
    newfile.write(comp_names + '\t' + species1 + '\t')
    newfile.write('\t'.join(list(map(lambda x: str(x), species_data[comp_names][0]))) + '\n')
    # write comp_names, species_2 and corresponding sites    
    newfile.write(comp_names + '\t' + species2 + '\t')
    newfile.write('\t'.join(list(map(lambda x: str(x), species_data[comp_names][1]))) + '\n')
    
# close file after writing
newfile.close()    
    
