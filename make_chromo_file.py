# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 14:21:28 2016

@author: RJovelin
"""

# use this script to make a file of chromosome names matched to scaffold/contig names

import os
import sys
import numpy as np
from scipy import stats
import math


# get option to use stingent or inclusive CNV definitions
CNVFilter = sys.argv[1]
assert CNVFilter in ['stringent', 'inclusive'], 'should use appropriate option'


# make a list of fasta files
files = [i for i in os.listdir() if i[-3:] == '.fa']
print(len(files))

# make a dictionary of scaffold: chromosome
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
            assert scaffold not in chromos, 'scaffold is already recorded'            
            chromos[scaffold] = LG
    infile.close()
print('done matching chromosome names', len(chromos))
                

# get CNVR coordinates {CNVR: [chromo, start, end]}
if CNVFilter == 'stringent':
    CNVFile = 'Stringent.Gain+Loss.hg19.2015-02-03.txt'
elif CNVFilter == 'inclusive':
    CNVFile = 'Inclusive.Gain+Loss.hg19.2015-02-03.txt'

CNVCoord = {}
infile = open(CNVFile)   
# skip header, read file
infile.readline()
for line in infile:
    line = line.rstrip()
    if line != '':
        line = line.split()
        # get chromo, start, end, state, region
        chromo, start, end, state, CNVR = line[0], int(line[1]) -1, int(line[2]), line[3], line[4]
        # check that state is CNV
        if state == 'CNV':
            CNVCoord[CNVR] = [chromo, start, end]
infile.close()
print('got CNVR coordinates', len(CNVCoord))

      
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
print('matched gene ID to gene name', len(GeneIDToGeneName))

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
print('matched mRNA ID to gene ID', len(mRNAToGene))

# make a dict with gene ID as key and a list of corresponding mRNA ID {gene ID: [rna ID]}
GeneTomRNA = {}
for rna in mRNAToGene:
    if mRNAToGene[rna] in GeneTomRNA:
        GeneTomRNA[mRNAToGene[rna]].append(rna)
    else:
        GeneTomRNA[mRNAToGene[rna]] = [rna]
print('matched Gene ID with all their mRNA ID', len(GeneTomRNA))

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
            LG = line[0]
            # get chromo name (eg chr1)
            if LG in chromos:
                chromo = chromos[LG]
                # get orientation
                orientation = line[6]
                # get start, end positions 0-based
                start = int(line[3]) -1
                end = int(line[4])
                mRNACoord[rna_id] = [chromo, start, end, orientation]  
# close file
infile.close()
print('extracted mRNA coordinates', len(mRNACoord))    
    
    
# record all overlaps between mRNAs and CNVR
overlap = {}
# count the number of cnv and non-cnv mRNAs
a, b, c = 0, 0, len(mRNACoord)
# record the CNV status of all mRNAs {rna ID: CNV status}
mRNACNV = {}
# find all mRNAs affected by CNVR
for rna in mRNACoord:
    c -= 1
    # get mRNA coord
    rna_chromo, rna_start, rna_end = mRNACoord[rna][0], mRNACoord[rna][1], mRNACoord[rna][2]
    # set boolean
    FoundCNV = False    
    # loop through CNVR
    for CNVR in CNVCoord:
        # get chromo, start and end
        cnv_chromo, cnv_start, cnv_end  = CNVCoord[CNVR][0], CNVCoord[CNVR][1], CNVCoord[CNVR][2]
        if cnv_chromo == rna_chromo:
            overlapping = len(set(range(cnv_start, cnv_end)).intersection(set(range(rna_start, rna_end))))
            if overlapping != 0:
                # record overlap
                if rna in overlap:
                    overlap[rna].append(overlapping)
                else:
                    overlap[rna] = [overlapping]
                # update boolean
                FoundCNV = True
                # exit loop
                break
    # check if the mrna overlaps with a CNV region
    if FoundCNV == True:
        # rna is found in CNVR
        mRNACNV[rna] = 'CNV'
        a += 1
    else:
        # rna not found in any of the CNVR
        mRNACNV[rna] = 'not_CNV'
        b += 1
    print('cnv: {0}, non-cnv: {1}, remaining: {2}'.format(a, b, c), sep = '\t', end = '\r')   
             

size = []
for rna in overlap:
    size.extend(overlap[rna])
if len(size) != 0:
    print(min(size), max(size), np.median(size), np.mean(size))