# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 14:10:59 2016

@author: RJovelin
"""


# use this script to compare the number target sites per nucleotide between CNV and non-CNV genes
# for single human populations 

# use Agg backend on server without X server
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rc
rc('mathtext', default='regular')
# import modules
import numpy as np
from scipy import stats
import math
import os
import sys
# import custom modules
from CNV_miRNAs import *


# usage PlotTargetsSingleHumanPop.py [options]
# [3UTR/5UTR/CDS] choose the region to analyse
# -[targetscan/miranda] : the predictor algorithm used to predict target sites


# get the region to consider to predict target sites [3UTR or 5UTr or CDS]
domain = sys.argv[1]
print(domain)
# get predictor from command [targetscan or miranda]
predictor = sys.argv[2]
print(predictor)

# keep genes on assembled chromsomes
keep_valid_chromos = True
chromos = 'valid_chromos'
# use all CNVs
cnv_length = 'CNV_all_length'

# get the CNV files with single human population references
CNV_file = 'GRCh37_hg19_variants_2016-05-15.txt'

# get the file with 3'UTR length
UTR_file = 'H_sapiens_3UTR_length_' + chromos + '.txt'
print(UTR_file)
# get targetscan sequence input file
targetscan_seq_input_file = 'H_sapiens_' + domain + '_' + chromos + '_targetscan.txt'
print(targetscan_seq_input_file)
# get the outputfile with predicted target sites
predicted_target_file = 'H_sapiens_' + domain + '_' + chromos + '_predicted_sites_' + predictor + '.txt'    
print(predicted_target_file)
# make a dictionary with {gene :[targets, seq_length, normalized_targets]}
if predictor == 'targetscan':
    predicted_targets = parse_targetscan_output(targetscan_seq_input_file, predicted_target_file, 'all')
elif predictor == 'miranda':
    predicted_targets = parse_miranda_output(targetscan_seq_input_file, predicted_target_file, 'all')

# make a dictionary {reference: pubmedid} for studies reported in the CGV
references = get_DGV_references(CNV_file)

# make a list of pubmed IDs for single population studies
pubmed = [['Suktitipat_et_al_2014', '25118596'], ['John_et_al_2014', '26484159'], ['Thareja_et_al_2015', '25765185'], ['Alsmadi_et_al_2014', '24896259']]

# check that reference names correspond to the pubmed IDs
infile = open(CNV_file)
for line in infile:
    for i in range(len(pubmed)):
        if pubmed[i][0] in line:
            assert pubmed[i][1] in line, 'study ref does not correspond to expected pubmed ID'
infile.close()


# create a set of CNV genes for each study {study: {gene1, gene2...}}
StudiesCNV = {}
for i in range(len(pubmed)):
    # get the set of CNV genes for that study
    CNV_genes = get_human_CNV_genes_single_study(CNV_file, pubmed[i][0], 'all')
    print(pubmed[i][0], len(CNV_genes))
    # populate dict
    StudiesCNV[pubmed[i][0]] = CNV_genes
print('extracted CNV genes for each study')

# get synonym names for all genes {gene name : [list of synonyms]}
synonyms = get_synonyms('H_sapiens.gff3')
print('got gene synonyms', len(synonyms))    
# get the CDS sequences of the longest mRNAs for each gene {gene : sequence}
CDS_seq = extract_CDS_sequences('H_sapiens.gff3', 'H_sapiens_genome.txt', 'H_sapiens_valid_chromos.txt', keep_valid_chromos)
print('extracted CDS sequences', len(CDS_seq))  


# create a dict {study: {gene: CNV status}}    
CNV_status = {}
for study in StudiesCNV:
    # initialize dict
    CNV_status[study] = {}
    # loop over genes in CDS seq
    for gene in CDS_seq:
        # set boolean
        is_cnv = False
        # ask if gene in CNV genes
        if gene in StudiesCNV[study] or gene.upper() in StudiesCNV[study]:
            # gene is CNV, add gene and status to inner dict
            CNV_status[study][gene] = 'CNV'
        else:
            # ask if any of the gene synonyms are in CNV genes
            for name in synonyms[gene]:
                # check if gene in CNV
                if name in StudiesCNV[study] or name.upper() in StudiesCNV[study]:
                    # update boolean variable
                    is_cnv = True
            # check if gene in CNV
            if is_cnv == True:
                CNV_status[study][gene] = 'CNV'
            elif is_cnv == False:
                CNV_status[study][gene] == 'not_CNV'
print('sorted genes according to CNV status')    
    
    
# make a dictionary {study: {gene: [targets, seq_length, normalized_targets, CNV_status]}}
CNVTargets = {}
for study in CNV_status:
    # initialize inner dict
    CNVTargets[study] = {}
    # loop over genes with CNV status for given study
    for gene in CNV_status[study]:
        # populate with list of targets
        CNVTargets[study][gene] = list(predicted_targets[gene])
        # add CNV status
        CNVTargets[study][gene].append(CNV_status[study][gene])
print('combined targets and CNV status')

    
# make a dict with the number of CNV genes for each study
CNVNum = {}
for study in CNV_status:
    CNVNum[study] = 0
    for gene in CNV_status[study]:
        if CNV_status[study][gene] == 'CNV':
            CNVNum[study] += 1
print('got CNV gene counts for each study')
for study in CNVNum:
    print(study, CNVNum[study])


##### continue here


            # parse the summary table into a list
            regulation = compare_miRNA_regulation('Temp_summary_targets.txt')
    
            # write regulation to file
            newfile.write(study + '\t')
            newfile.write('\t'.join(list(map(Gstr, regulation))) + '\n')
    
            print('done writing regulation for {0}'.format(study))

