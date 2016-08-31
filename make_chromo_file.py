# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 14:21:28 2016

@author: RJovelin
"""

# use this script to make a file of chromosome names matched to scaffold/contig names

import os
import numpy as np
from scipy import stats
import math


# make a list of fasta files
files = [i for i in os.listdir() if i[-3:] == '.fa']

# make a dictionary of chromosome: scaffold
chromos = {}

# loop over fasta files
for filename in files:
    infile = open(filename)
    # look for sequence headers
    for line in infile:
        if line.startswith('>'):
            line = line.rstrip().split('|')
            scaffold = line[3]
            LG = line[-1]
            LG = LG[:LG.index(',')]
            LG = LG.replace('Homo sapiens', '')
            LG = LG.replace(' ', "")
            LG = LG.replace('omosome', '')
            if LG in chromos:
                chromos[LG].append(scaffold)
            else:
                chromos[LG] = [scaffold]
    infile.close()
print('done matching chromosome names', len(chromos))
                

# get CNVR coordinates for Stringent and Inclusive filters
StringentCoord = {}
infile = open('Stringent.Gain+Loss.hg19.2015-02-03.txt')   
# skip header, read file
infile.readline()
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split()
        # get chromo, start, end, state
        chromo, start, end, state = line[0], int(line[1]) -1, int(line[2]), line[3]
        # check that state is CNV
        if state == 'CNV':
            StringentCoord = [chromo, start, end]
infile.close()
print('got Stringent CNVR coordinates', len(StringentCoord))

InclusiveCoord = {}
infile = open('Inclusive.Gain+Loss.hg19.2015-02-03.txt')   
# skip header, read file
infile.readline()
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split()
        # get chromo, start, end, state
        chromo, start, end, state = line[0], int(line[1]) -1, int(line[2]), line[3]
        # check that state is CNV
        if state == 'CNV':
            InclusiveCoord = [chromo, start, end]
infile.close()
print('got Inclusive CNVR coordinates', len(InclusiveCoord))
       
# get the gene coordinates, keeping the longest mRNA per gene
GFF_file = 'ref_GRCh37.p5_top_level.gff3'       
       
# get a dict with transcript ID : gene ID pair
TranscriptToGeneID = {}
infile = open(GFF_file)
# loop over file
for line in infile:
    # get line with mRNA
    if 'mRNA' in line:
        line = line.rstrip().split('\t')
        if line[2] == 'mRNA':
            # get chromo
            chromo = line[0]
            # parse descriptor string
            description = line[-1].split(';')
            # loop over string in description list
            for i in range(len(description)):
                # find and extract gene ID
                if 'GeneID:' in description[i]:
                    gene = description[i]
                # find and extract transcript id
                if 'transcript_id=' in description[i]:
                    TS = description[i]
            # further parse transcript ID and gene ID
            if ',' in gene:
                # check if comma happens before or after geneID
                if gene.count(',') == 1 and gene.index(',') < gene.index('GeneID:'):
                    gene = gene[gene.index('GeneID:') + 7:]
                else:
                    gene = gene[gene.index('GeneID:') + 7: gene.index(',', gene.index('GeneID:'))]                    
            else:
                gene = gene[gene.index('GeneID:') + 7: ]
            if ',' in TS:
                TS = TS[TS.index('transcript_id=') + 14: TS.index(',', TS.index('transcript_id'))]
            else:
                TS = TS[TS.index('transcript_id=') + 14: ]
            TranscriptToGeneID[TS] = gene
# close file
infile.close()
print('matched transcripts to parent genes', len(TranscriptToGeneID))


# make a dict to match genes to transcripts {gene_id: [list of transcripts]}
GeneToTranscript = {}
for TS in TranscriptToGeneID:
    # check if gene is key in dict
    if TranscriptToGeneID[TS] in GeneToTranscript:
        # add transcript to list
        GeneToTranscript[TranscriptToGeneID[TS]].append(TS)
    else:
        # create a list of transcript
        GeneToTranscript[TranscriptToGeneID[TS]] = [TS]
print('matched genes to transcripts', len(GeneToTranscript))

# match gene ID to gene name  {gene ID : gene name}
GeneIDToGeneName = {}
# open file for reading
infile = open(GFF_file)
# loop over file
for line in infile:
    # get line with mRNA
    if 'gene' in line:
        line = line.rstrip().split('\t')
        if line[2] == 'gene':
            # get chromo
            chromo = line[0]
            # parse descriptor string
            # separate on ';' because line doesn't always have the same structure
            description = line[-1].split(';')
            # loop over strings in list, find and extract GeneID
            for i in range(len(description)):
                if 'GeneID:' in description[i]:
                    ID = description[i]
                if 'Name=' in description[i]:
                    name = description[i]
            # parse ID and name
            if ',' in ID:
                # check if comma happens before or after geneID
                if ID.count(',') == 1 and ID.index(',') < ID.index('GeneID:'):
                    ID = ID[ID.index('GeneID:') + 7:]
                else:
                    ID = ID[ID.index('GeneID:') + 7: ID.index(',', ID.index('GeneID:'))]
            else:
                ID = ID[ID.index('GeneID:') + 7:]
            if ',' in name:
                name = name[name.index('Name=') + 5: name.index(',', name.index('Name='))]
            else:
                name = name[name.index('Name=') + 5:]
            GeneIDToGeneName[ID] = name
# close file
infile.close()


# match RNA ID to gene ID {RNA ID : gene ID} 
mRNAToGene = {}
# open file for reading
infile = open(GFF_file, 'r')
# loop over file
for line in infile:
    # find lines with mRNA
    if 'mRNA' in line:
        line = line.rstrip().split('\t')
        if line[2] == 'mRNA':
            # get chromo
            chromo = line[0]
            # extract rna id
            rna_id = line[-1][line[-1].index('ID=') + 3: line[-1].index(';')]
            # parse description line and extract gene ID
            description = line[-1].split(';')
            for i in range(len(description)):
                # find and extract gene ID
                if 'GeneID:' in description[i]:
                    gene = description[i]
            # further parse transcript ID and gene ID
            if ',' in gene:
                # check if comma happens before or after geneID
                if gene.count(',') == 1 and gene.index(',') < gene.index('GeneID:'):
                    gene = gene[gene.index('GeneID:') + 7:]
                else:
                    gene = gene[gene.index('GeneID:') + 7: gene.index(',', gene.index('GeneID:'))]                    
            else:
                gene = gene[gene.index('GeneID:') + 7: ]                    
            # populate dict 
            mRNAToGene[rna_id] = gene
# close file
infile.close()
    



# use this function to get the mRNA ID of each CDS ID
def match_CDS_to_mRNA(GFF_file, valid_chromo_file, keep_valid_chromos):
    '''
    (file, file, bool) -> dict
    Return a dictionary with CDS ID : mRNA ID pairs for all chromosomes
    (including MT, unplaced and unlocalized) if keep_valid_chromos 
    is False or for assembled nuclear chromosomes if True
    '''
    
    # make a set of valid chromos
    valid_chromos = get_valid_chromos(valid_chromo_file)
    
    # create dict {CDS ID : mRNA ID}
    CDS = {}
    
    # open file for reading
    infile = open(GFF_file, 'r')
    # loop over file
    for line in infile:
        # find lines with CDS
        if 'CDS' in line:
            line = line.rstrip().split('\t')
            # check that line corresponds to CDS info
            if line[2] == 'CDS':
                # get chromo
                chromo = line[0]
                # extract CDS ID
                cds_id = line[-1][line[-1].index('ID=') + 3: line[-1].index(';')]
                # extract parent mRNA ID
                rna_id = line[-1][line[-1].index('Parent=') + 7: line[-1].index(';', line[-1].index('Parent='))]
                # populate dict (each cds id maps to a single rna id)
                # ignore genes that are not linked to a mRNA (they map to gene segment) 
                if 'rna' in rna_id:
                    # check if need to consider only valid chromos
                    if keep_valid_chromos == True:
                        # check that gene in valid chromo
                        if chromo in valid_chromos:
                            CDS[cds_id] = rna_id
                    elif keep_valid_chromos == False:
                        # record all genes, regardless of linkage
                        CDS[cds_id] = rna_id                
                
    # close file
    infile.close()
    return CDS


# use this function to map mRNA ID to their CDS id
def match_mRNA_to_CDS(GFF_file, valid_chromo_file, keep_valid_chromos):
    '''
    (file, file, bool) -> dict
    Return a dictionary with mRMA ID : CDS ID pairs for all chromosomes
    (including MT, unplaced and unlocalized) if keep_valid_chromos 
    is False or for assembled nuclear chromosomes if True
    '''
    
    # get a dict with CDS ID : mRNA ID pairs
    CDS = match_CDS_to_mRNA(GFF_file, valid_chromo_file, keep_valid_chromos)
            
    # reverse dict
    mRNAs = {}
    for cds_id in CDS:
        mRNAs[CDS[cds_id]] = cds_id
        
    return mRNAs


# get the coordinates of all mRNAs
# create a dict {rna_id: [chromo, start, end , orientation]}
mRNACoord = {}
# open file for reading
infile = open(GFF_file, 'r')
# loop over file
for line in infile:
    # ignore lines that do not contain mRNA
    if 'mRNA' in line:
        line = line.rstrip().split('\t')
        # check that line correspond to a a mRNA
        if line[2] == 'mRNA':
            # extract RNA ID
            rna_id = line[-1][line[-1].index('ID=') + 3: line[-1].index(';')]
            # get chromo
            chromo = line[0]
            # get orientation
            orientation = line[6]
            # get start, end positions 0-based
            start = int(line[3]) -1
            end = int(line[4])
            mRNACoord[rna_id] = [chromo, start, end, orientation]  
# close file
infile.close()
    

 
    
    # Note that a few rnas correspond to gene_ids that map to the same gene name
    # keep a single rna corresponding to a single gene id that nap to a single gene name
    # Rnor: 49; Hsa: 0; Ptr: 11; Cfa: 18; Bt: 94; Mmul: 335; Gga: 34; Mmus: 2

    # create a set to store the longest mRNA for each gene (ie, a single mRNA per gene)
    longest = set()
    
    # get the mRNA coordinates {mRNA : [chromo, start, end, orientation]}
    mRNA_coord = get_mRNA_coord(GFF_file, valid_chromo_file, keep_valid_chromos)
    # get the {mRNA ID : gene ID} pairs
    mRNA_gene = match_mRNA_to_gene(GFF_file, valid_chromo_file, keep_valid_chromos)
    # get the {gene ID: gene name} pairs
    gene_names = geneID_to_gene_name(GFF_file, valid_chromo_file, keep_valid_chromos)
    # get the {rna ID : CDS ID} pairs
    mRNA_CDS = match_mRNA_to_CDS(GFF_file, valid_chromo_file, keep_valid_chromos)
       
    # compute the length of each mRNA
    mRNA_length = {}
    for rna in mRNA_coord:
        L = mRNA_coord[rna][2] - mRNA_coord[rna][1]
        mRNA_length[rna] = L
    # create a dictionary with gene name as key : list of [L, rna]
    genes = {}
    for rna in mRNA_length:
        # ignore rna_id not mapping to CDS ID
        if rna in mRNA_CDS:
            # get gene id
            gene_id = mRNA_gene[rna]
            # get gene name
            name = gene_names[gene_id]
            # populate dict
            if name in genes:
                # add [mRNA_length, mRNA]
                genes[name].append([mRNA_length[rna], rna])
            else:
                genes[name] = [[mRNA_length[rna], rna]]
            
    # sort mRNA according to their length
    for name in genes:
        genes[name].sort()
        # get the longest mRNA
        rna = genes[name][-1][-1]
        longest.add(rna)
    
    return longest
        
        
    
    longest = get_longest_mRNA(GFF_file, valid_chromo_file, keep_valid_chromos)
        
    # get the mRNA coordinates {rna_id: [chromo, start, end , orientation]}
    mRNA_coord = get_mRNA_coord(GFF_file, valid_chromo_file, keep_valid_chromos)
    
    # get the CDS coordinates {cds_id: [chromo, orientation, [[start, end], [start, end]]}
    CDS_coord = get_CDS_coord(GFF_file, valid_chromo_file, keep_valid_chromos)
    
    # get the cds_id of each rna_id {rna_id : cds id}
    rna_to_cds = match_mRNA_to_CDS(GFF_file, valid_chromo_file, keep_valid_chromos)
        
    # get the gene id of each rna id {RNA ID : gene ID}
    rna_to_gene = match_mRNA_to_gene(GFF_file, valid_chromo_file, keep_valid_chromos)
    
    # get the gene name for each gene id {gene ID : gene name pairs}
    gene_names = geneID_to_gene_name(GFF_file, valid_chromo_file, keep_valid_chromos)
        





