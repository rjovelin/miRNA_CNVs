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
print('matched gene to transcripts', len(TranscriptToGeneID))


# make a dict to match genes to transcripts {gene_id: [list of transcripts]}
GeneToTranscript = {}


def (GFF_file, valid_chromo_file, keep_valid_chromos):
    '''
    (file, file, bool) -> dict
    Return a dictionary with gene ID : list of of transcript pairs for all
    chromosomes (including MT, unplaced and unlocalized) if keep_valid_chromos 
    is False or for assembled nuclear chromosomes if True
    '''
    
    # create a dict of transcript ID : gene ID
    transcripts = transcript_to_geneID(GFF_file, valid_chromo_file, keep_valid_chromos)
        
    # reverse dictionary
    genes = {}
    for TS in transcripts:
        # check if gene is key in dict
        if transcripts[TS] in genes:
            # add transcript to list
            genes[transcripts[TS]].append(TS)
        else:
            # create a list of transcript
            genes[transcripts[TS]] = [TS]
            
    return genes


# use this function to get a dict with gene ID : gene name pair
def geneID_to_gene_name(GFF_file, valid_chromo_file, keep_valid_chromos):
    '''
    (file, file, bool) -> dict
    Return a dictionary with gene ID : gene name pairs for all
    chromosomes (including MT, unplaced and unlocalized) if keep_valid_chromos 
    is False or for assembled nuclear chromosomes if True
    '''
    
    # create a set of valid chromosomes
    valid_chromos = get_valid_chromos(valid_chromo_file)
        
    # create a dict
    genes = {}    
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
                
                # check if need to consider only valid chromos
                if keep_valid_chromos == True:
                    # check that gene in valid chromo
                    if chromo in valid_chromos:
                        genes[ID] = name
                elif keep_valid_chromos == False:
                    # record all genes, regardless of linkage
                    genes[ID] = name
                
    # close file
    infile.close()
    
    return genes



# use this function to match a RNA ID with a gene ID
def match_mRNA_to_gene(GFF_file, valid_chromo_file, keep_valid_chromos):
    '''
    (file, file, bool) -> dict
    Return a dictionary with RNA ID : gene ID pairs for all
    chromosomes (including MT, unplaced and unlocalized) if keep_valid_chromos 
    is False or for assembled nuclear chromosomes if True
    '''
    
    # make a set of valid chromos
    valid_chromos = get_valid_chromos(valid_chromo_file)
       
    # create a dict {RNA ID : gene ID}
    rnas = {}
    
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
                 
                # check if need to consider only valid chromos
                if keep_valid_chromos == True:
                    # check that gene in valid chromo
                    if chromo in valid_chromos:
                        rnas[rna_id] = gene
                elif keep_valid_chromos == False:
                    # record all genes, regardless of linkage
                    rnas[rna_id] = gene

    # close file
    infile.close()
    
    return rnas


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


# use this function to get the coordinates of each mRNA
def get_mRNA_coord(GFF_file, valid_chromo_file, keep_valid_chromos):
    '''
    (file, file, bool) -> dict
    Return a dictionary with RNA ID and a list with coordinates for all
    chromosomes (including MT, unplaced and unlocalized) if keep_valid_chromos 
    is False or for assembled nuclear chromosomes if True
    '''
    
    # make a set of valid chromos
    valid_chromos = get_valid_chromos(valid_chromo_file)
    
    # create a dict {rna_id: [chromo, start, end , orientation]}
    rnas = {}
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
                
                # check if need to consider only valid chromos
                if keep_valid_chromos == True:
                    # check that gene in valid chromo
                    if chromo in valid_chromos:
                        rnas[rna_id] = [chromo, start, end, orientation]    
                elif keep_valid_chromos == False:
                    # record all genes, regardless of linkage
                    rnas[rna_id] = [chromo, start, end, orientation]  
                
    # close file
    infile.close()
    
    return rnas

 
# use this function to get the CDS coordinates of each mRNA
def get_CDS_coord(GFF_file, valid_chromo_file, keep_valid_chromos):
    '''
    (file, file, bool) -> dict
    Return a dictionary with CDS ID as key and a list of CDS coordinates as value
    for all chromosomes (including MT, unplaced and unlocalized) if 
    keep_valid_chromos is False or for assembled nuclear chromosomes if True
    '''
    
    # make a set of valid chromos
    valid_chromos = get_valid_chromos(valid_chromo_file)

    # create a dict to store CDS coordinates
    # {cds_ID: [chromo, orientation, [[start, end], [start, end]]}
    CDS_coord = {}
    
    # open file for reading
    infile = open(GFF_file, 'r')
    
    # loop over file
    for line in infile:
        # create a boolean to record CDS coordinates only under some condition
        record_CDS = False
        # ignore line that do not contain CDS
        if 'CDS' in line:
            line = line.rstrip().split('\t')
            if line[2] == 'CDS':
                # get chromo
                chromo = line[0]
                # get orientation
                orientation = line[6]
                # get start, end position 0-based
                start = int(line[3]) - 1
                end = int(line[4])
                # extract the cds id
                cds_id = line[-1][line[-1].index('ID=') + 3: line[-1].index(';')]
                # ignore CDS that are not linked to a mRNA (these are linked to gene segment)
                parent = line[-1][line[-1].index('Parent=') + 7: line[-1].index(';', line[-1].index('Parent='))]
                if 'rna' in parent:
                    # check if need to consider only valid chromos
                    if keep_valid_chromos == True:
                        # check that gene in valid chromo
                        if chromo in valid_chromos:
                            record_CDS = True
                    elif keep_valid_chromos == False:
                        # record all genes, regardless of linkage
                        record_CDS = True                   
                    # check if CDS can be recorded (ie, in valid chromo or all chromo)
                    if record_CDS == True:
                        # check if cds_id is key in dict
                        if cds_id in CDS_coord:
                            # add CDS coordinates to list
                            CDS_coord[cds_id][2].append([start, end])
                        else:
                            # create list with chromo, orientation and [CDS start, CDS end]
                            CDS_coord[cds_id] = [chromo, orientation, [[start, end]]]
    
    # close file after reading
    infile.close()
        
    # make sure that coordinates are sorted according to their order on chromo
    for cds_id in CDS_coord:
        CDS_coord[cds_id][2].sort()
    
    return CDS_coord



# use this function to get the rna_id of the longest mRNA for each gene
def get_longest_mRNA(GFF_file, valid_chromo_file, keep_valid_chromos):
    '''
    (file, file, bool) -> set
    Return a set with the longest mRNA for each gene (ie. a single mRNA per gene)
    for all chromosomes (including MT, unplaced and unlocalized) if 
    keep_valid_chromos is False or for assembled nuclear chromosomes if True
    '''
    
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
        
        
    
# use this function to get the 5' UTR coordinates of the longest mRNA of each gene
def get_5UTR_coord(GFF_file, valid_chromo_file, keep_valid_chromos):
    '''
    (file, file, bool) -> dict
    Return a dictionary with gene as key and a list of 5' UTR coordinates
    for all chromosomes (including MT, unplaced and unlocalized) if 
    keep_valid_chromos is False or for assembled nuclear chromosomes if True
    '''
    
    # create a dict with UTR_coord {gene_name: [chromo, start, end, orientation]}
    UTR_coord = {}
    
    # keep the longest mRNA of each gene
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
        
    # loop over mRNA in mRNA cood
    for rna in mRNA_coord:
        # check if rna is among the longest
        if rna in longest:
            # get mRNA orientation 
            orientation = mRNA_coord[rna][-1]
            # get mRNA chromo
            chromo = mRNA_coord[rna][0]
            # check orientation
            if orientation == '+':
                # 5'UTR start = mRNA start
                UTR_start = mRNA_coord[rna][1] 
                # 5'UTR end = start position of 1st CDS
                UTR_end = CDS_coord[rna_to_cds[rna]][2][0][0]
            elif orientation == '-':
                # UTR start is end of last CDS
                UTR_start = CDS_coord[rna_to_cds[rna]][2][-1][1]
                # UTR end is end of mRNA
                UTR_end = mRNA_coord[rna][2]
            
            # check that chromo and orientation match between mRNA and CDS
            assert CDS_coord[rna_to_cds[rna]][0] == chromo, 'chromo doesn\'t match between mRNA and CDS'
            assert CDS_coord[rna_to_cds[rna]][1] == orientation, 'orientation doesn\'t match between mRNA and CDS'
            # get corresponding gene id
            gene_id = rna_to_gene[rna]
            # get gene name
            gene = gene_names[gene_id]
            # check that gene name is unique
            assert gene not in UTR_coord, 'gene name already recorded'
            # populate dict with chromo, UTR_start, UTR_end, orientation
            UTR_coord[gene] = [chromo, UTR_start, UTR_end, orientation]
    
    return UTR_coord



# use this function to get the 3' UTR coordinates of the longest mRNA of each gene
def get_3UTR_coord(GFF_file, valid_chromo_file, keep_valid_chromos):
    '''
    (file, file, bool) -> dict
    Return a dictionary with gene as key and a list of 3' UTR coordinates
    for all chromosomes (including MT, unplaced and unlocalized) if 
    keep_valid_chromos is False or for assembled nuclear chromosomes if True
    '''
    
    # create a dict with UTR_coord {gene_name: [chromo, start, end, orientation]}
    UTR_coord = {}
    
    # keep the longest mRNA of each gene
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
        
    # loop over mRNA in mRNA cood
    for rna in mRNA_coord:
        # check if rna is among the longest
        if rna in longest:
            # get mRNA orientation 
            orientation = mRNA_coord[rna][-1]
            # get mRNA chromo
            chromo = mRNA_coord[rna][0]
            # check orientation
            if orientation == '+':
                # 3'UTR start is end of last CDS
                UTR_start = CDS_coord[rna_to_cds[rna]][2][-1][1]
                # 3'UTR end = end of mRNA
                UTR_end = mRNA_coord[rna][2]
            elif orientation == '-':
                # UTR_start is start of mRNA
                UTR_start = mRNA_coord[rna][1]
                # UTR_end is start of 1st CDS
                UTR_end = CDS_coord[rna_to_cds[rna]][2][0][0]
                
            # check that chromo and orientation match between mRNA and CDS
            assert CDS_coord[rna_to_cds[rna]][0] == chromo, 'chromo doesn\'t match between mRNA and CDS'
            assert CDS_coord[rna_to_cds[rna]][1] == orientation, 'orientation doesn\'t match between mRNA and CDS'
            # get corresponding gene id
            gene_id = rna_to_gene[rna]
            # get gene name
            gene = gene_names[gene_id]
            # check that gene name is unique
            assert gene not in UTR_coord, 'gene name already recorded'
            # populate dict with chromo, UTR_start, UTR_end, orientation
            UTR_coord[gene] = [chromo, UTR_start, UTR_end, orientation]
    
    
    return UTR_coord


# use this function to extract the CDS sequence of the longest mRNAs
def extract_CDS_sequences(GFF_file, genome_fasta, valid_chromo_file, keep_valid_chromos):
    '''
    (file, file, file, bool) -> dict
    Return a dictionary with gene name : coding sequence pairs for gene on all
    chromosomes (including MT, unplaced and unlocalized) if keep_valid_chromos
    is False or for assembled nuclear chromosomes if True.
    '''
    
    # create a dictionary to store the CDS sequences {gene name : CDS sequence}
    CDS_sequences = {}
    
    # convert genome file to dict
    genome = convert_genome_fasta(genome_fasta)
    
    # get the coordinates of all CDS {cds_ID: [chromo, orientation, [[start, end], [start, end]]}
    CDS_coord = get_CDS_coord(GFF_file, valid_chromo_file, keep_valid_chromos)
    
    # need to keep the CDS coordinates of the longest mRNAs 
    longest = get_longest_mRNA(GFF_file, valid_chromo_file, keep_valid_chromos)
        
    # get the rna_id of each cds id
    cds_to_rna = match_CDS_to_mRNA(GFF_file, valid_chromo_file, keep_valid_chromos)
    
    # get the gene id of each rna id {RNA ID : gene ID}
    rna_to_gene = match_mRNA_to_gene(GFF_file, valid_chromo_file, keep_valid_chromos)
    
    # get the gene name for each gene id {gene ID : gene name pairs}
    gene_names = geneID_to_gene_name(GFF_file, valid_chromo_file, keep_valid_chromos)
    
    # loop of cds_id:
    for cds_id in CDS_coord:
        # get rna_id
        rna_id = cds_to_rna[cds_id]
        # check if rna_id in longest
        if rna_id in longest:
            # get gene if
            gene_id = rna_to_gene[rna_id]
            # get gene name
            gene = gene_names[gene_id]
            # get chromo
            chromo = CDS_coord[cds_id][0]
            # get orientation
            orientation = CDS_coord[cds_id][1]
            # create a sequence string
            sequence = ''
            # loop over CDS coordinates, grab exon sequence, update CDS sequence
            for pair in CDS_coord[cds_id][2]:
                # get start and end positions 0-based
                start = pair[0]
                end = pair[1]
                # get exon sequence
                exon = genome[chromo][start:end].upper()
                sequence += exon
            # take the reverse complement if orientation is -
            if orientation == '-':
                sequence = reverse_complement(sequence)
            # update dict with gene name : sequence pair
            CDS_sequences[gene] = sequence
            
    return CDS_sequences        


# use this function to extract the UTR sequences of the longest mRNA for each gene
def extract_UTR_sequences(GFF_file, genome_fasta, valid_chromo_file, UTR_type, keep_valid_chromos):
    '''
    (file, file, file, str, bool) -> dict
    Return a dictionary with gene name : UTR sequence pairs for 5'UTR or 3'UTR
    depending on UTR_type and for genes on all chromosomes (including MT,
    unplaced and unlocalized) if keep_valid_chromos is False or for assembled
    nuclear chromosomes if True.
    '''

    # get UTR coordinates {gene_name: [chromo, start, end, orientation]}
    # check UTR type
    if UTR_type == '5UTR':
        # get 5'UTR coordinates
        UTR_coord = get_5UTR_coord(GFF_file, valid_chromo_file, keep_valid_chromos)
    elif UTR_type == '3UTR':
        UTR_coord = get_3UTR_coord(GFF_file, valid_chromo_file, keep_valid_chromos)
    
    # convert genome to dict
    genome = convert_genome_fasta(genome_fasta)
    
    # create a dict {gene name : UTR sequence}
    UTR_sequences = {}

    # loop over gene name in UTR coord
    for gene in UTR_coord:
        # get chromo
        chromo = UTR_coord[gene][0]
        # get start and end position 0-based
        start = UTR_coord[gene][1]
        end = UTR_coord[gene][2]
        # get orientation
        orientation = UTR_coord[gene][-1]
        # grab UTR sequence
        sequence = genome[chromo][start:end].upper()
        # take reverse complement is orienation is -
        if orientation == '-':
            sequence = reverse_complement(sequence)
        # populate dict
        UTR_sequences[gene] = sequence
        
    return UTR_sequences
       
       
       
       
       
       