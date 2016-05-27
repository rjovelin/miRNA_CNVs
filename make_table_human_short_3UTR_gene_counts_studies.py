# -*- coding: utf-8 -*-
"""
Created on Sun Oct 11 11:31:09 2015

@author: Richard
"""

# usage make_table_human_short_3UTR_gene_counts_studies.py [True/False] [long_CNVs/all_CNVs] DGV_file


# determine the number of human genes with short (< 7b) 3'UTR in each study for a given version of the DGV

from CNV_miRNAs import *
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
print(keep_valid_chromos, chromos)

# get the option to call a CNV if CNV length > 1 Kb (long_CNVs)
# or to include all CNVs regardless of length (all_CNVs)
long_CNV = sys.argv[2]
if long_CNV == 'all_CNVs':
    cnv_length = 'CNV_all_length'
    CNV_size = 'all'
elif long_CNV == 'long_CNVs':
    cnv_length = 'CNV_greater_1Kb'
    CNV_size = 'long'
print(long_CNV, cnv_length, CNV_size)

# get DGV file from the command
CNV_file = sys.argv[3]
print(CNV_file)

# note: only CNVs on valid chromos are reported in DGV, so if all chromos are
# used it may introduce a bias by calling non CNV genes genes that cannot be detected

# get UTR_file
UTR_file = 'H_sapiens_3UTR_length_' + chromos + '.txt'
print(UTR_file)

# get release version
if '2013' in CNV_file:
    release_version = 'GRCh37_' + CNV_file[CNV_file.rindex('_') + 1: -7]
else:
    release_version = 'GRCh37_' + CNV_file[CNV_file.rindex('_') + 1: -10]
print(release_version)    

# get outputfile
outputfile = 'Human_Counts_Short3UTR_' + cnv_length + '_' + chromos + '_' + release_version + '.txt'
print(outputfile)

# open file for writing
newfile = open(outputfile, 'w')

# write headers to file
newfile.write('# Number of Human genes with short (< 7 bp) 3\'UTR\n')
newfile.write('\t'.join(['Study', 'PubMedID', 'Total', 'CNV', 'non-CNV']) + '\n')

# get the dictionaries of {reference: pubmedid} for studies reported in the CGV
references = get_DGV_references(CNV_file)

# get synonym names for all genes {gene name : [list of synonyms]}
synonyms = get_synonyms('H_sapiens.gff3')
    
# loop over study
for study in references:
    print(study)
    # get the set of CNV genes corresponding to that study
    CNV_genes = get_human_CNV_genes_single_study(CNV_file, study, CNV_size)
    
    print('# CNV genes', len(CNV_genes))    
    
    # sort genes based on 3' UTR length
    UTR_length = sort_genes_3UTR_length(UTR_file)
    print(len(UTR_length))

    # get the CNV status of all short #' UTR genes
    # create a dict {gene: CNV status}    
    CNV_status = {}
    
    # loop over gene in UTR_length
    for gene in UTR_length:
        # set boolean
        is_cnv = False
        # check UTR length
        if UTR_length[gene] == 'short':
            # ask if gene in CNV genes
            if gene in CNV_genes or gene.upper() in CNV_genes:
                # gene is CNV, add gene and status to dict
                CNV_status[gene] = 'CNV'
            else:
                # ask if any of the gene synonyms are in CNV genes
                for name in synonyms[gene]:
                    # check if in CNV genes
                    if name in CNV_genes or name.upper() in CNV_genes:
                        # update boolean variable
                        is_cnv = True
                # check if gene is CNV
                if is_cnv == True:
                    CNV_status[gene] = 'CNV'
                elif is_cnv == False:
                    CNV_status[gene] = 'not_CNV'
    
    print('# short genes with CNV status', len(CNV_status))
    
    # count total number of genes with short 3'UTR
    total_short = 0
    # count CNV genes with short UTR
    cnv_short = 0
    # count non-CNV genes with short UTR
    non_cnv_short = 0    
    
    # loop over genes in UTR_length
    for gene in UTR_length:
        if UTR_length[gene] == 'short':
            total_short += 1
            # check CNV status
            if CNV_status[gene] == 'CNV':
                cnv_short += 1
            elif CNV_status[gene] == 'not_CNV':
                non_cnv_short += 1
            
    print('total short', total_short)
    print('CNV short', cnv_short)
    print('non_CNV', non_cnv_short)
    
    assert total_short == cnv_short + non_cnv_short, 'sum cnv and non-cnv short is not equal to total short'
        
    # write results to file
    newfile.write('\t'.join([study, references[study], str(total_short), str(cnv_short), str(non_cnv_short)]) + '\n')
    
    
# close file after writing
newfile.close()

