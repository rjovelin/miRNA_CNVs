# -*- coding: utf-8 -*-
"""
Created on Sat Sep 26 10:31:09 2015

@author: Richard
"""

import numpy as np
from scipy import stats
import math



# use this function to get the 
def read_genome(filename):
    '''
    (file) -> dict
    Take a file containing a single chromosome sequence in Fasta format and return
    a dictionary with the sequence header as key and a single string sequence as value
    Precondition: the file contains a single sequence (a single header with >)
    '''
    
    # initiate dict
    genome = {}    
    # open file for reading
    infile = open(filename, 'r')
    # get the sequence header
    header = infile.readline().rstrip()
    # get rid of the > sign
    header = header[1:]
    # get the entire sequence
    sequence = infile.read().rstrip()
    # close file
    infile.close()
    # get read of end of line
    sequence = sequence.replace('\n', '')
    # populate dict
    genome[header] = sequence
    # close file
    infile.close()

    return genome    
    

# use this function to translate a DNA sequence
def cds_translate(cds):
    '''
    (str) -> str
    Translate a coding sequence into a protein sequence according to the standard genetic code

    >>> cds_translate('ATGGCCATGGCGCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA')
    MAMAPRTEINSTRING*
    >>> cds_translate('ATGTACTAA')
    MY*
    '''

    genetic_code = {'TTT': 'F', 'CTT': 'L', 'ATT': 'I', 'GTT': 'V',
                   'TTC': 'F', 'CTC': 'L', 'ATC': 'I', 'GTC': 'V',
                   'TTA': 'L', 'CTA': 'L', 'ATA': 'I', 'GTA': 'V',
                   'TTG': 'L', 'CTG': 'L', 'ATG': 'M', 'GTG': 'V',
                   'TCT': 'S', 'CCT': 'P', 'ACT': 'T', 'GCT': 'A',
                   'TCC': 'S', 'CCC': 'P', 'ACC': 'T', 'GCC': 'A',
                   'TCA': 'S', 'CCA': 'P', 'ACA': 'T', 'GCA': 'A',
                   'TCG': 'S', 'CCG': 'P', 'ACG': 'T', 'GCG': 'A',
                   'TAT': 'Y', 'CAT': 'H', 'AAT': 'N', 'GAT': 'D',
                   'TAC': 'Y', 'CAC': 'H', 'AAC': 'N', 'GAC': 'D',
                   'TAA': '*', 'CAA': 'Q', 'AAA': 'K', 'GAA': 'E',
                   'TAG': '*', 'CAG': 'Q', 'AAG': 'K', 'GAG': 'E',
                   'TGT': 'C', 'CGT': 'R', 'AGT': 'S', 'GGT': 'G',
                   'TGC': 'C', 'CGC': 'R', 'AGC': 'S', 'GGC': 'G',
                   'TGA': '*', 'CGA': 'R', 'AGA': 'R', 'GGA': 'G',
                   'TGG': 'W', 'CGG': 'R', 'AGG': 'R', 'GGG': 'G'}
   

    CDS = cds.upper()
    protein = ''

    for i in range(0, len(CDS), 3):
        codon = CDS[i:i+3]
        if codon not in genetic_code:
            protein += 'X'
        else:
            protein += genetic_code[codon]

    return protein

# use this function to convert a fasta file into a dictionary
def convert_fasta(fasta):
    '''
    (file) -> dict
    Take a file with fasta sequences and return a dictionnary with
    sequence ID as key and single string sequence as value
    '''
    # convert nematode genomes into a dictionnary
    genome = {}
    infile = open(fasta, 'r')
    for line in infile:
        line = line.rstrip()
        if line != '':
            if line.startswith('>'):
                genome[line[1:]] = ""
                seq_name = line[1:]
            else:
                genome[seq_name] += line
    infile.close
    return genome


# use this function to create a dict of chromo : sequence pairs
def convert_genome_fasta(genome_fasta):
    '''
    (file) -> dict
    Take the genome fasta file with sequence ID in GenBank format and return
    a dictionary of chromo : sequence pairs    
    '''
    
    # convert to the fasta file to a dict
    genome_gbk = convert_fasta(genome_fasta)
    
    # create a dict chromo : sequence pairs
    genome = {}
    # loop over genome, extract chromo, populate new dict
    for seq_id in genome_gbk:
        # parse the sequence id
        chromo = seq_id.split('|')
        chromo = chromo[3]
        genome[chromo] = genome_gbk[seq_id]
        
    return genome
        

# use this function to take the reverse complement of a DNA sequence
def reverse_complement(dna):
    '''
    (str) -> (str)
    Return the reverse complementary sequence of string dna

    >>> reverse_complement('atcg')
    'cgat'
    '''

    valid_bases = {'A', 'T', 'C', 'G'}

    dna2 = dna.upper()
    dna_comp = ''
    for i in dna2:
        if i == 'A':
            dna_comp += 'T'
        elif i == 'T':
            dna_comp += 'A'
        elif i == 'C':
            dna_comp += 'G'
        elif i == 'G':
            dna_comp += 'C'
        elif i not in valid_bases:
            dna_comp += 'N'

    reverse_comp_dna = ''
    for i in reversed(dna_comp):
        reverse_comp_dna += i

    if dna.islower():
        reverse_comp_dna = reverse_comp_dna.lower()
        
    return reverse_comp_dna



# use this function to make a set of valid chromosomes
def get_valid_chromos(valid_chromo_file):
    '''
    (file) -> set
    Return a set of valid chromosomes corresponding to the nuclear genome
    '''
    
    # Note: The DataBase of Genomic Variants removes any variants (and genes)
    # not mapping to the autosomes and sex chromsomes
    # keep only autosomes and sex chromosomes for all species and ignore unplaced,
    # unlocalized, mitochondrial genes
    
    # create a set of chromos
    valid = set()
    
    # open file for reading
    infile = open(valid_chromo_file, 'r')
    # loop over file
    for line in infile:
        line = line.rstrip()
        if line != '':
            valid.add(line)
    # close file after reading
    infile.close()
            
    return valid



# use this function to get a dict with transcript ID : gene ID pair
def transcript_to_geneID(GFF_file, valid_chromo_file, keep_valid_chromos):
    '''
    (file, file, bool) -> dict
    Return a dictionary with transcript ID : gene ID pairs for all
    chromosomes (including MT, unplaced and unlocalized) if keep_valid_chromos 
    is False or for assembled nuclear chromosomes if True
    '''
    
    # create a set of valid chromosomes
    valid_chromos = get_valid_chromos(valid_chromo_file)
    
    # create a dict
    transcripts = {}    
    # open file for reading
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
                      
                # check if need to consider only valid chromos
                if keep_valid_chromos == True:
                    # check that gene in valid chromo
                    if chromo in valid_chromos:
                        transcripts[TS] = gene
                elif keep_valid_chromos == False:
                    # record all genes, regardless of linkage
                    transcripts[TS] = gene
    # close file
    infile.close()
    
    return transcripts



# use this function to make a dict of gene : list of transcripts
def geneID_to_transcript(GFF_file, valid_chromo_file, keep_valid_chromos):
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
        
        

# use this function to group miRNAs by families
def mirna_family(mature_fasta):
    '''
    (file) -> dict
    Take a fasta file of mature miRNA sequences and return a dictionary with
    seed as key and a list of miRNA sharing the same seed, ie from the same family
    '''
    
    # convert fasta file to dict
    mirnas = convert_fasta(mature_fasta)
    # reverse the dictionary
    families = {}
    # loop over mirnas
    for mature in mirnas:
        # grab seed
        seed = mirnas[mature][1:8]
        # check if seed in dict, add mature to list
        if seed in families:
            families[seed].append(mature)
        else:
            families[seed] = [mature]

    return families            


# use this function to make a dict of mirna name : seed
def mirna_to_seed(mature_fasta):
    '''
    (file) -> dict
    Take a fasta file of mature miRNA sequences and return a dictionary with
    mirna name : seed pairs
    '''
    
    # make a dict of seed : list of mirnas
    seeds = mirna_family(mature_fasta)
    
    # create a dict of mirna {name : seed}
    mirnas = {}
    
    # reverse dictionary of seed : list of name, remove miRBAse ID in name
    for seq in seeds:
        # loop over the names
        for name in seeds[seq]:
            # remove miRBAse ID, remove species code
            name = name[name.index('-') + 1:name.index(' ')]
            # populate dict
            mirnas[name] = seq
            
    return mirnas




# use this function to get the the mirnas families shared by 2 species
def find_conserved_mirna_families(mature_fasta1, mature_fasta2):
    '''
    (file, file) -> set
    Take the fasta files of mature sequences for species 1 and for species 2
    and return a set of mirna seeds shared by the 2 species
    '''
    
    # get the mirna families in each species
    fam_sp1 = mirna_family(mature_fasta1)
    fam_sp2 = mirna_family(mature_fasta2)
    
    # create sets of seeds
    seeds_sp1 = set(i for i in fam_sp1)
    seeds_sp2 = set(i for i in fam_sp2)
    
    return seeds_sp1.intersection(seeds_sp2)


# use this function to get the miRNA families conserved among a group of species
def find_shared_mirna_families(mature_fasta_files):
    '''
    (list) -> set
    Take a list of fasta files of mature mirna sequences for a group of species
    and return a set with all the seeds found in each species (ie. only the seeds
    that are shared by all species)
    '''
    
    # create a set of shared families
    shared = set()    
    
    # initialize set by adding seeds of first species
    fam = mirna_family(mature_fasta_files[0])
    shared = set(k for k in fam)    
    
    # remove from shared any seed that is not found in the other species
    for i in range(1, len(mature_fasta_files)):
        # get the seeds of the species
        fam = fam = mirna_family(mature_fasta_files[i])
        seeds = set(k for k in fam)
        # remove seeds not present in species
        to_remove = set()
        to_remove = shared - seeds
        for j in to_remove:
            shared.remove(j)
    
    return shared
    

# use this function to get all the mirna families found in a group of species
def seeds_species_group(mature_fasta_files):
    '''
    (list) -> set
    Take a list of fasta files of mature mirna sequences for a group of species
    and return a set with all the seeds found in the group (ie. found at least once
    in any given species)
    '''
    
    # make a set of all seeds present in the group
    total = set()
    # loop over list of files
    for i in range(len(mature_fasta_files)):
        # get the mirna families in each species
        fam = mirna_family(mature_fasta_files[i])
        # create a set of seeds for that species
        seeds = set(j for j in fam)
        # add each seed to total seeds
        for j in seeds:
            total.add(j)
            
    return total
    
    
   
# use this function to find the mirna families unique to a given species
def find_species_specific_mirna_families(mature_fasta1, mature_fasta_files):
    '''
    (file, list) -> set
    Take the fasta file of mature mirna sequences in species 1, and a list
    of fasta files of mature mirna sequences for a group of species and return
    a set of seeds that are unique to species 1 (ie. all seeds in present in
    species and not found in any of the other species)
    '''
    
    # get the seeds of mirnas in species 1
    fam_sp1 = mirna_family(mature_fasta1)
    # create set of seeds
    seeds_sp1 = set(i for i in fam_sp1)
    
    # create a set of seeds present in the group of species
    # check if focal species in the group
    if mature_fasta1 in mature_fasta_files:
        # create a list of fasta files without the focal species
        group = mature_fasta_files[:]
        group.remove(mature_fasta1)
        print(mature_fasta1, len(mature_fasta_files), len(group))
    else:
        group = mature_fasta_files[:]
        
    # get all the seeds present in the group of species
    seeds_group = seeds_species_group(group)
    print(len(seeds_sp1), len(seeds_group))
    
    return seeds_sp1 - seeds_group
    
    

# use this function to generate a set a set of all human gene names
def get_all_human_genes(GFF_file):
    '''
    (file) -> set
    Return a set with all human gene names
    '''
    
    # create a set
    genes = set()    
    # open file for reading
    infile = open(GFF_file)
    # loop over file
    for line in infile:
        # get line with gene
        if 'gene' in line:
            line = line.rstrip().split('\t')
            if line[2] == 'gene':
                # parse descriptor string
                # separate on ';' because line doesn't always have the same structure
                description = line[-1].split(';')
                # loop over strings in list, find and extract GeneID
                for i in range(len(description)):
                    if 'Name=' in description[i]:
                        name = description[i]
                # name
                if ',' in name:
                    name = name[name.index('Name=') + 5: name.index(',', name.index('Name='))]
                else:
                    name = name[name.index('Name=') + 5:]
                # add name to set
                genes.add(name)
                                
    # close file
    infile.close()
    
    return genes


# use this function to get a set of human CNV genes
def get_human_CNV_genes(CNV_file, CNV_size):
    '''
    (file, str) -> set
    Returns a set of human CNV genes, for CNV of any length if CNV_size is "all",
    or for CNV of minimal length of 1 KB if CNV_size is "long"
    '''

    # make a set to store the CNV genes
    # Note that the affected genes are sometimes line[-1] when no sample ID is not provided
    # or line[-2] when sample ID is provided: take all and get CNV genes status by asking for membership

    # create a set of CNV genes
    CNV_genes = set()

    # open file for reading
    infile = open(CNV_file, 'r')
    # skip header
    infile.readline()
    # loop over file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split('\t')
            # select only CNV variants
            if line[4] == 'CNV':
                # get start and position 0-based
                start = int(line[2]) - 1
                end = int(line[3])
                # get CNV length
                CNV_L = end - start
                if ',' in line[-2]:
                    genes = line[-2].split(',')
                    for gene in genes:
                        # check if size limit
                        if CNV_size == 'long' and CNV_L >= 1000:
                            # record only CNV greater than 1 KB
                            CNV_genes.add(gene)
                        elif CNV_size == 'all':
                            # record all CNV genes
                            CNV_genes.add(gene)
                elif ',' not in line[-2]:
                    # check if size limit
                    if CNV_size == 'long' and CNV_L >= 1000:
                        # record only CNV greater than 1 KB
                        CNV_genes.add(line[-2])
                    elif CNV_size == 'all':
                        # record all CNV genes
                        CNV_genes.add(line[-2])
                  
                if ',' in line[-1]:
                    genes = line[-1].split(',')
                    for gene in genes:
                        # check if size limit
                        if CNV_size == 'long' and CNV_L >= 1000:
                            # record only CNV greater than 1 KB
                            CNV_genes.add(gene)
                        elif CNV_size == 'all':
                            # record all CNV genes
                            CNV_genes.add(gene)
                elif ',' not in line[-1]:
                    if CNV_size == 'long' and CNV_L >= 1000:
                        # record only CNV greater than 1 KB
                        CNV_genes.add(line[-1])
                    elif CNV_size == 'all':
                        # record all CNV genes
                        CNV_genes.add(line[-1])
                    
    infile.close()
    return CNV_genes



# use this function to get references supporting CNV genes in DGV 
def get_DGV_references(CNV_file):
    '''
    (file) -> dict
    Return a dictionary with study name : PubMed_ID pairs for all
    studies reported in the DataBase of Genomic Variants
    '''
    
    # create a dict {reference: pubmedid}
    references = {}
    
    # open file for reading
    infile = open(CNV_file, 'r')
    # skip header
    infile.readline()
    # loop over file
    for line in infile:
        # record only studies reporting CNV variants
        if 'CNV' in line:
            line = line.rstrip().split('\t')
            # check variant type
            if line[4] == 'CNV':
                references[line[6]] = line[7]
                
    # close file
    infile.close()
    
    return references
   


# use this function to extract CNV genes for a given study
def get_human_CNV_genes_single_study(CNV_file, study, CNV_size):
    '''
    (file, str, str) -> set
    Returns a set of CNV genes for a given study in the DataBase of Genomic
    Variants for all CNVs if CNV_size is "all", or for CNV of minimal length
    of 1 KB if CNV_size is "long"
    '''

    # make a set to store the CNV genes
    # Note that the affected genes are sometimes line[-1] when no sample ID is not provided
    # or line[-2] when sample ID is provided: take all and filter out the non valid genes

    infile = open(CNV_file, 'r')
    CNV_genes = set()
    infile.readline()
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split('\t')
            if line[4] == 'CNV' and line[6] == study:
                # get CNV position
                start = int(line[2]) - 1
                end = int(line[3])
                # get CNV length
                CNV_L = end - start
                if ',' in line[-2]:
                    genes = line[-2].split(',')
                    for gene in genes:
                        # check size limit
                        if CNV_size  == 'long' and CNV_L > 1000:
                            # record only CNV greater than 1 Kb
                            CNV_genes.add(gene)
                        elif CNV_size == 'all':
                            # record all CNVs
                            CNV_genes.add(gene)
                elif ',' not in line[-2]:
                    # check size limit
                    if CNV_size  == 'long' and CNV_L > 1000:
                        # record only CNV greater than 1 Kb
                        CNV_genes.add(line[-2])
                    elif CNV_size == 'all':
                        # record all CNVs
                        CNV_genes.add(line[-2])
                    
                if ',' in line[-1]:
                    genes = line[-1].split(',')
                    for gene in genes:
                        # check size linit
                        if CNV_size  == 'long' and CNV_L > 1000:
                            # record only CNV greater than 1 Kb
                            CNV_genes.add(gene)
                        elif CNV_size == 'all':
                            # record all CNVs
                            CNV_genes.add(gene)
                elif ',' not in line[-1]:
                    # check size limit
                    if CNV_size  == 'long' and CNV_L > 1000:
                        # record only CNV greater than 1 Kb
                        CNV_genes.add(line[-1])
                    elif CNV_size == 'all':
                        # record all CNVs
                        CNV_genes.add(line[-1])


    infile.close()
    return CNV_genes



# use this function to generate the targetscan UTR input file
def make_targetscan_seq_input_file(GFF_file, genome_fasta, valid_chromo_file, keep_valid_chromos, gene_region, species, outputfile):
    
    '''
    (file, file, file, bool, str, str, file) -> file
    Save the sequences in the input format of TargetScan for "5UTR", "3UTR" or
    "CDS" depending on the type of gene_region and for genes on all chromosomes
    (including MT, unplaced and unlocalized) if keep_valid_chromos is False
    or for assembled nuclear chromosomes if True.
    '''
        
    # make a dict with species : NCBI species ID
    species_ID = {'Hsa': '9606', 'Ptr': '9598', 'Mmul': '9544', 'Mmus': '10090',
                  'Rno': '10116', 'Bta': '9913', 'Cfa': '9615', 'Gga': '9031'}
    
    # make a dict of gene name : UTR sequence pairs
    # check gene region to consider
    if gene_region == '5UTR':
        sequences = extract_UTR_sequences(GFF_file, genome_fasta, valid_chromo_file, '5UTR', keep_valid_chromos)
    elif gene_region == '3UTR':
        sequences = extract_UTR_sequences(GFF_file, genome_fasta, valid_chromo_file, '3UTR', keep_valid_chromos)        
    elif gene_region == 'CDS':
        sequences = extract_CDS_sequences(GFF_file, genome_fasta, valid_chromo_file, keep_valid_chromos)
        
    # open file for writing
    newfile = open(outputfile, 'w')
    
    # loop over genes in sequences
    for gene in sequences:
        # check if sequence if empty of not
        if sequences[gene] != '':
            # format for targetscan input is tab delimited: gene, species_id, mRNA sequence
            newfile.write(gene + '\t' + species_ID[species] + '\t' + sequences[gene].upper().replace('T', 'U') + '\n')
        
    # close file after writing
    newfile.close()  

   
# use this function to slice the 7bp seed from the mature sequence 
def grab_seed(mature_seq):
    '''
    (str) -> str
    Slice the miRNA mature sequence to return the 7bp seed motif
    '''
    seed = mature_seq[1:8]
    return seed


# use this function to generate a dict of seeds and list of mirnas pairs
def seed_mirnas(mature_fasta):
    '''
    (file) -> dict
    Take a fasta file of miRNA mature sequences and return a dictionnary with
    seed as key and a list of mirnas from the same family (sharing the same seed)
    as value
    '''
    
    # convert fasta file to 
    mirnas = convert_fasta(mature_fasta)
    
    # create a dict of seed : [mir1, mir2]
    seeds = {}
    
    # loop over mature sequences
    for name in mirnas:
        # get seed sequence
        seed_seq = grab_seed(mirnas[name])
        # clean up the mirna name
        mirna_name = name[4:name.index(' ')]
        # populate dict
        if seed_seq in seeds:
            seeds[seed_seq].append(mirna_name)
        else:
            seeds[seed_seq] = [mirna_name]
            
    return seeds
    
    
# use this function to generate the TargetScan miRNA family input file
def make_targetscan_mirfam_input_file(mature_fasta, species, outputfile):
    '''
    (file) -> file
    Take a fasta file of miRNA mature sequences and generate the miRNA family
    input file for TargetScan    
    '''
    
    # make a dict with species : NCBI species ID
    species_ID = {'Hsa': '9606', 'Ptr': '9598', 'Mmul': '9544', 'Mmus': '10090',
                  'Rno': '10116', 'Bta': '9913', 'Cfa': '9615', 'Gga': '9031'}
    
    # create a dict of seed sequence and list of mirna pairs
    seeds = seed_mirnas(mature_fasta)
    
    # open file for writing
    newfile = open(outputfile, 'w')
    
    # loop over the seed sequences
    for motif in seeds:
        # record only 1 mirna per family
        # grab first mirna
        # write mirna, seed, species ID to file
        newfile.write(seeds[motif][0] + '\t' + motif + '\t' + species_ID[species] + '\n')
    
    # close file after writing
    newfile.close()



# use this function to sort genes according to their 3'UTR length
def sort_genes_3UTR_length(gene_UTR_length_file, L):
    '''
    (file, int) -> dict
    Return a dictionary with gene name : 3'UTR length description pairs, 
    short if length < L bp and long if length > L bp    
    '''

    # create a dict {gene : "short" (or "long")}
    genes = {}
    
    # open file for reading
    infile = open(gene_UTR_length_file, 'r')
    # skip header
    infile.readline()
    # loop over file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            # populate dict
            if int(line[1]) < L:
                genes[line[0]] = 'short'
            else:
                genes[line[0]] = 'long'
    # close file
    infile.close()
    
    return genes


# use this function to sort genes according to their CNV status
def sort_genes_CNV_status(gene_CNV_status_file):
    '''
    (file) -> dict
    Return a dictionary with gene : CNV status pairs
    Note: CNV status is already called based on CNV length, and 
    '''
    
    # create dict to {gene : CNV status}
    genes ={}
    # open file for reading
    infile = open(gene_CNV_status_file, 'r')
    # skip header
    infile.readline()
    # loop over file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            genes[line[0]] = line[1]
    # close file
    infile.close()
    
    return genes
    

# use this function to make a set of CNV genes
def retrieve_CNV_genes(CNV_file):
    '''
    (file) -> set
    Take the file with CNV genes obtained from the published data and a return 
    a set if CNV genes
    Precondition: CNV genes have already been called based on the length of CNVs
    with indication of CNV length in file name (all_length, or greater than 1 Kb)
    '''
    
    # create a set of CNV genes
    CNV_genes = set()
    
    # open file for reading
    infile = open(CNV_file, 'r')
    # loop over file
    for line in infile:
        line = line.rstrip()
        if line != '':
            # add gene to set
            CNV_genes.add(line)
            
    # close file
    infile.close()
    
    return CNV_genes


# use this function to sort genes according to their CNV status
def get_genes_CNV_status(GFF_file, genome_fasta, valid_chromo_file, keep_valid_chromos, CNV_file):
    '''
    (file, file, file, bool, file) -> dict
    Return a dictionary with gene name : CNV status (CNV or not_CNV) for gene
    on all chromosomes (including MT, unplaced and unlocalized) if keep_valid_chromos
    is False or for assembled nuclear chromosomes if True
    '''
    
    # create a dict {gene: CNV status}    
    CNV_status = {}
    
    # get the CDS sequences of the longest mRNAs for each gene {gene : sequence}
    CDS_seq = extract_CDS_sequences(GFF_file, genome_fasta, valid_chromo_file, keep_valid_chromos)

    # get synonym names for all genes {gene name : [list of synonyms]}
    synonyms = get_synonyms(GFF_file)
    
    # get the set of published CNV genes
    CNV_genes = retrieve_CNV_genes(CNV_file)
    
    # loop over gene in CDS seq
    for gene in CDS_seq:
        # set boolean
        is_cnv = False
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
    
    return CNV_status


# use this function to find synonymous names for all genes
def get_synonyms(GFF_file):
    '''
    (file) -> dict
    Return a dictionary with gene name : list of synonymous names pairs
    '''
    
    # create a dict {gene name : [list of synonyms]}
    gene_synonyms = {}    
    
    # open file for reading
    infile = open(GFF_file)
    # loop over file
    for line in infile:
        # get line with mRNA
        if 'gene' in line:
            line = line.rstrip().split('\t')
            if line[2] == 'gene':
                # separate on ';' because line doesn't always have the same structure
                description = line[-1].split(';')
                # find gene name
                for i in range(len(description)):
                    if 'Name=' in description[i]:
                        name = description[i]
                # parse name
                if ',' in name:
                    name = name[name.index('Name=') + 5: name.index(',', name.index('Name='))]
                else:
                    name = name[name.index('Name=') + 5:]
                # set up synonym variable
                synonym = ''
                # find synonyms
                for i in range(len(description)):
                    if 'gene_synonym' in description[i]:
                        synonym = description[i]
                # check if gene has synonyms
                if synonym == '':
                    # no synonym, add gene name to itself to capture all genes
                    gene_synonyms[name] = [name]
                else:
                    # parse synonyms
                    synonym = synonym[synonym.index('gene_synonym=') + 13:]
                    # check if synonym contains extra text
                    if ';' in synonym:
                        synonym = synonym[:synonym.index(';')]
                    # get all synonymous names
                    if ',' in synonym:
                        # parse on ',' to get all names
                        synonym = synonym.split(',')
                        # populate dict with gene : list pair
                        gene_synonyms[name] = synonym
                    else:
                        # populate dict with the single synonym
                        gene_synonyms[name] = [synonym]
                                
    # close file
    infile.close()
    
    return gene_synonyms

                
  
# use this function to get the length of the region used to predict target sites
def get_domain_length_from_targetscan_input(targetscan_seq_input_file):
    '''
    (file) -> dict
    Return a dictionary with gene : sequence length pairs
    '''
    
    # create dict {gene : sequence length}
    genes = {}
    
    # open file for reading
    infile = open(targetscan_seq_input_file, 'r')
    # loop over file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split('\t')
            # populate dict
            genes[line[0]] = len(line[2])
            
    # close file
    infile.close()
    
    return genes
    


      
# use this function to create dicts from the targetscan output
def parse_targetscan_output(targetscan_seq_input_file, predicted_targets, miRNAs, *seeds_mature):
    '''
    (file, file, str, *set, *file) -> dict
    Take the targetscan input sequence file, the targetscan output, 
    and return a dictionary with gene as key and a list with the number of
    predicted target sites, sequence length and number of target sites
    normalized by sequence length for all miRNAs or for conserved miRNAs only.
    If miRNAs = conserved, the a set of conserved seeds and the fasta file with 
    mature miRNAs is obtained from the optional arguments
    Precondition: files are generated to include valid or all chromos
    '''
    
    # get seed set from optional parameter tuple
    if miRNAs == 'conserved':
        seeds = seeds_mature[0]
        # make a dict of mirna {name : seed} pairs
        mature_fasta = seeds_mature[1]
        mir_names = mirna_to_seed(mature_fasta)
    
    # get the length of the sequences used to predict target sites {gene : seq_length}
    genes_length = get_domain_length_from_targetscan_input(targetscan_seq_input_file)
    
    # create a dict to store the target sites {gene: set(target1, target2, target3)} 
    target_counts = {}
    
    # record the number of target sites from the targetscan output    
    
    # open file for reading
    infile = open(predicted_targets, 'r')
    # skip header
    infile.readline()
    # go through file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split('\t')
            # get gene name
            gene = line[0]
            # get site defined as "start:end:site_type"
            # some sites may be bound by different miRNA families
            start = line[3]
            end = line[4]
            site_type = line[8]
            site = start + ':' + end + ':' + site_type
            # get mirna binding to site
            mirna = line[1]
            # check if consider all or conserved mirna families
            if miRNAs == 'conserved':
                # check that mirna seed in set of seeds
                if mir_names[mirna] in seeds:
                    # populate dict
                    if gene not in target_counts:
                        # initialize value with empty set
                        target_counts[gene] = set()
                    # add target site to set
                    target_counts[gene].add(site)
            elif miRNAs == 'all':
                # record all sites
                # populate dict
                if gene not in target_counts:
                    # initialize value with empty set
                    target_counts[gene] = set()
                # add target site to set
                target_counts[gene].add(site)
                        
              
    # close file after reading
    infile.close()
    
    # create a dict to store the number of targets {gene: [N_sites, seq_length, N_sites/seq_length]}
    targets = {}
    
    # loop over genes with predicted targets in targetscan
    for gene in target_counts:
        # add number of target sites
        targets[gene] = [len(target_counts[gene])]
        # get the length of the sequence used to predict target sites
        seq_length = genes_length[gene]
        # sequence length should be > 7 bp for prediction
        assert seq_length >= 7, 'sites are predicted but sequence length is < 7 bp'
        # add sequence length to list
        targets[gene].append(seq_length)
        # add number of sites normalized by sequence length
        targets[gene].append(len(target_counts[gene]) / seq_length)
        
    # count 0 for genes that do not have any targets but that have a region > mininum length
    
    # loop over genes in targetscan seq input
    for gene in genes_length:
        # check that sites are not already recorded
        if gene not in targets:
            # targets were not predicted by targetscan
            # check that sequence length can lead to prediction
            if genes_length[gene] >= 7:
                # populate dict
                targets[gene] = [0, genes_length[gene], 0]
                
    return targets



# use this function to create dict from miranda outputs
def parse_miranda_output(targetscan_seq_input_file, predicted_targets, miRNAs, *seeds_mature):
    '''
    (file, file, str, *set, *file) -> dict
    Take the targetscan input sequence file, the miranda output, 
    and return a dictionary with gene as key and a list with the number of
    predicted target sites, sequence length and number of target sites
    normalized by sequence length for all miRNAs or for conserved miRNAs only.
    If miRNAs = conserved, the a set of conserved seeds and the fasta file with 
    mature miRNAs is obtained from the optional arguments
    Precondition: files are generated to include valid or all chromos    
    '''
    
    
    # get seed set from optional parameter tuple
    if miRNAs == 'conserved':
        seeds = seeds_mature[0]
        # make a dict of mirna {name : seed} pairs
        mature_fasta = seeds_mature[1]
        mir_names = mirna_to_seed(mature_fasta)    
    
    # get the length of the sequences used to predict target sites {gene : seq_length}
    genes_length = get_domain_length_from_targetscan_input(targetscan_seq_input_file)
        
    # create a dict to store the target sites {gene: set(target1, target2, target3)} 
    target_counts = {}
    
    # create a dict with {gene : sequence length} pairs
    target_length = {}    
    
    
    # open file for reading
    infile = open(predicted_targets, 'r')
    # go through file
    for line in infile:
        if line.startswith('>>'):
            line = line.rstrip().split('\t')
            # get mirna, get rid of '>>' sign
            mirna = line[0][2:]
            # get rid of species code in mirna name
            mirna = mirna[4:]
            # get target gene
            gene = line[1]
            # get sequence length
            seq_length = int(line[8])
            # populate dict with gene : sequence length pairs
            if gene not in target_length:
                target_length[gene] = seq_length
            # get positions
            positions = line[9].split()
            # check if consider all or conserved miRNA families
            if miRNAs == 'conserved':
                # check that mirna seed in set of seeds
                if mir_names[mirna] in seeds:
                    # record sites only if mirna seed in set of seeds
                    # populate dict
                    if gene not in target_counts:
                        # inititalise value with empty list
                        target_counts[gene] = set()
                    # add all sites to set
                    for pos in positions:
                        target_counts[gene].add(pos)
            elif miRNAs == 'all':
                # record all sites
                # populate dict
                if gene not in target_counts:
                    # initialize value with empty set
                    target_counts[gene] = set()
                for pos in positions:
                    target_counts[gene].add(pos)
                    
    # close file after reading
    infile.close()
    
    # create a dict to store the number of targets {gene: [N_sites, seq_length, N_sites/seq_length]}
    targets = {}
    
    # loop over genes with predicted targets in targetscan
    for gene in target_counts:
        # add number of target sites
        targets[gene] = [len(target_counts[gene])]
        # get the length of the sequence used to predict target sites
        seq_length = target_length[gene]
        # add sequence length to list
        targets[gene].append(seq_length)
        # add number of sites normalized by sequence length
        targets[gene].append(len(target_counts[gene]) / seq_length)
    
    
    # count 0 for genes that do not have any targets but that have a region > mininum length
    
    # loop over genes in targetscan seq input
    for gene in genes_length:
        # check that sites are not already recorded
        if gene not in targets:
            # targets were not predicted
            # check that sequence is greater than 6 bp
            if genes_length[gene] >= 7:
                # populate dict
                targets[gene] = [0, genes_length[gene], 0]

    return targets




###########################################

## use this function to create dicts from the targetscan output
#def filter_targetscan_output(targetscan_seq_input_file, predicted_targets, mature_fasta, seeds):
#    '''
#    (file, file, file, set) -> dict
#    Take the targetscan input sequence file, the targetscan output, the fasta file
#    with mature miRNAs, set of miRNA seeds and and return a dictionary with gene
#    as key and a list with the number of predicted target sites,
#    sequence length and number of target sites normalized by sequence length
#    only for mirnas for which the seed is included in the seed set
#    Precondition: files are generated to include valid or all chromos
#    '''
#    
#    # make a dict of mirna {name : seed} pairs    
#    mir_names = mirna_to_seed(mature_fasta)
#    
#    # get the length of the sequences used to predict target sites {gene : seq_length}
#    genes_length = get_domain_length_from_targetscan_input(targetscan_seq_input_file)
#    
#    # create a dict to store the target sites {gene: set(target1, target2, target3)} 
#    target_counts = {}
#    
#    # record the number of target sites from the targetscan output    
#    
#    # open file for reading
#    infile = open(predicted_targets, 'r')
#    # skip header
#    infile.readline()
#    # go through file
#    for line in infile:
#        line = line.rstrip()
#        if line != '':
#            line = line.split('\t')
#            # get gene name
#            gene = line[0]
#            # get mirna binding to site
#            mirna = line[1]
#            # record sites only if mirna seed in set if seeds
#            if mir_names[mirna] in seeds:
#                # get site defined as "start:end:site_type"
#                # some sites may be bound by different miRNA families
#                start = line[3]
#                end = line[4]
#                site_type = line[8]
#                site = start + ':' + end + ':' + site_type
#                # populate dict
#                if gene not in target_counts:
#                    # initialize value with empty set
#                    target_counts[gene] = set()
#                # add target site to set
#                target_counts[gene].add(site)
#            
#    # close file after reading
#    infile.close()
#    
#    # create a dict to store the number of targets {gene: [N_sites, seq_length, N_sites/seq_length]}
#    targets = {}
#    
#    # loop over genes with predicted targets in targetscan
#    for gene in target_counts:
#        # add number of target sites
#        targets[gene] = [len(target_counts[gene])]
#        # get the length of the sequence used to predict target sites
#        seq_length = genes_length[gene]
#        # sequence length should be > 7 bp for prediction
#        assert seq_length >= 7, 'sites are predicted but sequence length is < 7 bp'
#        # add sequence length to list
#        targets[gene].append(seq_length)
#        # add number of sites normalized by sequence length
#        targets[gene].append(len(target_counts[gene]) / seq_length)
#        
#    # count 0 for genes that do not have any targets but that have a region > mininum length
#    
#    # loop over genes in targetscan seq input
#    for gene in genes_length:
#        # check that sites are not already recorded
#        if gene not in targets:
#            # targets were not predicted by targetscan
#            # check that sequence length can lead to prediction
#            if genes_length[gene] >= 7:
#                # populate dict
#                targets[gene] = [0, genes_length[gene], 0]
#                
#    return targets    



###############################################################

           
# use this function to write summary table with target site counts and CNV status
def make_summary_table_target_sites(predicted_targets, gene_CNV_status_file, outputfile):
    '''
    (dict, file, file) -> file 
    Take the dictionary with target sites obtained by parsing the predicted
    targets outputfile, the file with CNV gene status, and write to outputfile
    a summary table with the number of targets, the length of the sequence used
    to predict targets, the number of targets normalized by sequence length,
    and the CNV status for each gene
    Precondition: files are generated to include valid or all chromos and
    CNVs > 1 Kb or all CNVs
    '''

    # get CNV_status
    CNV_status = sort_genes_CNV_status(gene_CNV_status_file)

    # open file for writing
    newfile = open(outputfile, 'w')
    # write header
    newfile.write('\t'.join(['Gene', 'N_targets', 'Sequence_length', 'N_targets_normalized', 'CNV_status']) + '\n')
    
    # loop over genes in targets
    for gene in predicted_targets:
        # write number of sites, sequence length and normalized number of targets to file
        newfile.write('\t'.join([gene, str(predicted_targets[gene][0]), str(predicted_targets[gene][1]), str(predicted_targets[gene][2])]) + '\t')
        # write CNV status
        newfile.write(CNV_status[gene] + '\n')
        
    # close file after writing
    newfile.close()



# use this function to compare the number of miRNA target sites between CMV and non-CNV genes
def compare_miRNA_regulation(summary_table):
    '''
    (file) -> list
    Take a file with the number of target sites, sequence length and normalized
    number of targets for each gene and return a list with average number of target sites, 
    average sequence and average normalized targets for CNV and non-CNV genes
    '''
    
    # make lists of sequence length for CNV and non CNV genes
    CNV_seq, nonCNV_seq = [], []

    # make lists of target sites for CNV and non CNV genes
    CNV_targets, nonCNV_targets = [], []

    # make lists of normalized targets for CNV and non CNV genes
    CNV_normalized, nonCNV_normalized = [], []    
    
    # open file for reading
    infile = open(summary_table, 'r')
    # skip header
    infile.readline()
    # loop over file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            # get number of targets
            N_targets = int(line[1])
            # get sequence length
            seq_length = int(line[2])
            # get normalized number of targets
            normalized = float(line[3])
            # get CNV status
            CNV_status = line[-1]
            # populate lists
            if CNV_status == 'CNV':
                CNV_seq.append(seq_length)
                CNV_targets.append(N_targets)
                CNV_normalized.append(normalized)
            elif CNV_status == 'not_CNV':
                nonCNV_seq.append(seq_length)
                nonCNV_targets.append(N_targets)
                nonCNV_normalized.append(normalized)
                
    # close file after reading
    infile.close()
    
    # compute mean sequence length, mean number of targets, mean number of normalized sites
    CNV_average_targets = np.mean(CNV_targets)
    nonCNV_average_targets = np.mean(nonCNV_targets)
    CNV_average_norm = np.mean(CNV_normalized)
    nonCNV_average_norm = np.mean(nonCNV_normalized)
    CNV_average_length = np.mean(CNV_seq)
    nonCNV_average_length = np.mean(nonCNV_seq)
    
    # compute standard error of the mean
    CNV_sem_targets = np.std(CNV_targets) / math.sqrt(len(CNV_targets))
    nonCNV_sem_targets = np.std(nonCNV_targets) / math.sqrt(len(nonCNV_targets))
    CNV_sem_norm = np.std(CNV_normalized) / math.sqrt(len(CNV_normalized))
    nonCNV_sem_norm = np.std(nonCNV_normalized) / math.sqrt(len(nonCNV_normalized))
    CNV_sem_length = np.std(CNV_seq) / math.sqrt(len(CNV_seq))
    nonCNV_sem_length = np.std(nonCNV_seq) / math.sqrt(len(nonCNV_seq))
    
    # get the P value of Wilcoxon sum rank tests between CNV and non CNV genes
    P_targets = stats.ranksums(CNV_targets, nonCNV_targets)[1]
    P_normalized = stats.ranksums(CNV_normalized, nonCNV_normalized)[1]
    P_length = stats.ranksums(CNV_seq, nonCNV_seq)[1]


    # compute correlation between number of target sites and sequence length
    # make a copy of the list of CNV targets
    all_targets = CNV_targets[:]
    # add all values of the targets of non-CNV genes
    all_targets.extend(nonCNV_targets)
    # copy sequence length of CNV targets
    all_seq = CNV_seq[:]
    # add all values of the length of non CNV genes
    all_seq.extend(nonCNV_seq)
    # get Spearman's correlation coefficient and P-value
    correl_coeff, P_correl = stats.spearmanr(all_targets, all_seq)
    
    return [len(CNV_targets), CNV_average_targets, CNV_sem_targets,
            len(nonCNV_targets), nonCNV_average_targets, nonCNV_sem_targets,
            P_targets, CNV_average_norm, CNV_sem_norm,
            nonCNV_average_norm, nonCNV_sem_norm, P_normalized, 
            CNV_average_length, CNV_sem_length, nonCNV_average_length,
            nonCNV_sem_length, P_length, correl_coeff, P_correl]



# use this function to parse the summary table of target sites
def parse_summary_table_targets(summary_table):
    '''
    (file) -> list of lists
    Take a file with the number of target sites, sequence length and normalized
    number of targets for each gene and return a list of lists with 
    number of target sites, sequence length, number of normalized targets
    for CNV and non-CNV genes
    '''
    
    # make lists of sequence length for CNV and non CNV genes
    CNV_seq, nonCNV_seq = [], []

    # make lists of target sites for CNV and non CNV genes
    CNV_targets, nonCNV_targets = [], []

    # make lists of normalized targets for CNV and non CNV genes
    CNV_normalized, nonCNV_normalized = [], []    
    
    # open file for reading
    infile = open(summary_table, 'r')
    # skip header
    infile.readline()
    # loop over file
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split()
            # get number of targets
            N_targets = int(line[1])
            # get sequence length
            seq_length = int(line[2])
            # get normalized number of targets
            normalized = float(line[3])
            # get CNV status
            CNV_status = line[-1]
            # populate lists
            if CNV_status == 'CNV':
                CNV_seq.append(seq_length)
                CNV_targets.append(N_targets)
                CNV_normalized.append(normalized)
            elif CNV_status == 'not_CNV':
                nonCNV_seq.append(seq_length)
                nonCNV_targets.append(N_targets)
                nonCNV_normalized.append(normalized)
                
    # close file after reading
    infile.close()
    
    
    return [CNV_targets, nonCNV_targets, CNV_seq, nonCNV_seq, CNV_normalized, nonCNV_normalized]    
    

# use this function to map miRNA accession number to mature accession numbers or mature names
def MatchmiRNAAccessionNumbers(MatureRecord, miRBaseFile = 'miRNA.dat'):
    '''
    (file, str) -> dict
    Take the file with miRNA information from miRBase and and return a dictionary
    with species name as key and a dictionary of mirna accession number : list of
    mature names or mature accessions pairs depending on the value the string
    variable MatureRecord    
    '''
        
    # create a dictionary {species: {mirna accession : [mature accession/names]}}
    miRNAs = {}    
    
    # loop over file with miRNA information
    infile = open(miRBaseFile)
    for line in infile:
        if line.startswith('AC'):
            # get the mirna accession number
            line = line.rstrip().split()
            miRNAAccession = line[1][:line[1].index(';')]
        elif line.startswith('DE'):
            # get the species name
            line = line.rstrip().split()
            species = line[1] + '_' + line[2]
            # initialize dict with species name
            if species not in miRNAs:
                miRNAs[species] = {}
            # initialize dict
            assert miRNAAccession not in miRNAs[species], 'mirna accession number already recorded'
            miRNAs[species][miRNAAccession] = []
        elif line.startswith('FT'):
            if MatureRecord == 'accession' and 'accession' in line:
                # record the mature accession number
                line = line.rstrip().split()
                mature = line[1][line[1].index('"')+1: -1]
                # add to list
                miRNAs[species][miRNAAccession].append(mature)
            elif MatureRecord == 'name' and 'product' in line:
                # record the mature name
                line = line.rstrip().split()
                mature = line[1][line[1].index('"')+1: -1]
                miRNAs[species][miRNAAccession].append(mature)
                
    infile.close()
    return miRNAs



# use this function to extract the expression of each miRNA
def miRBAsemiRNAExpression(ExpressionFile = 'mirna_read_count.txt'):
    '''
    (file) -> dict
    Take the miRBase file with expression level for each miRNA and return a
    dict with miRNA accession number : expression (in RPM) pairs    
    '''
    
    # create dict {accession: expression}
    miRNAExpression = {}
    # loop over file, match accession numbers with the corresponding expression level
    infile = open(ExpressionFile)
    for line in infile:
        line = line.rstrip()
        if line != '':
            line = line.split('\t')
            # get the accession number
            accession = line[1]
            # get the normalized expression in RPM
            expression = float(line[-1])
            assert accession not in miRNAExpression, 'mirna expression already recorded'
            miRNAExpression[accession] = expression
    infile.close()
    return miRNAExpression
    
    
# use this function to sort mirnas based on expression quartiles
def SortmiRNAQuartileExpression(species, miRNAExpression, miRNAs):
    '''
    (str, dict, dict) -> tuple
    Take a species name (Genus_species), the dict of mirna accession: expression
    level pairs and return a list of lists of mirna accessions sorted for the given
    species according to their expression
    '''
    
    # miRNAExpression is a dict in the form {accession: expression}
    # miRNAs is a dict with accession number paris for each species {species: {mirna accession : [mature accession/names]}}
    
    # create list of expression values
    expression_level = [miRNAExpression[mirna] for mirna in miRNAExpression if mirna in miRNAs[species]]
    # compute quartiles of expression values
    Q1 = np.percentile(expression_level, 25)
    Q2 = np.percentile(expression_level, 50)
    Q3 = np.percentile(expression_level, 75)
    
    # partition mirnas according to expression quartiles
    # create lists to store theta for different levels of expression
    lowExp, moderateExp, mediumExp, highExp = [], [], [], []
    # loop over mirna in theta dict
    for mirna in miRNAExpression:
        # check expression value
        if miRNAExpression[mirna] < Q1:
            lowExp.append(mirna)
        elif miRNAExpression[mirna] >= Q1 and miRNAExpression[mirna] < Q2:
            moderateExp.append(mirna)
        elif miRNAExpression[mirna] >= Q2 and miRNAExpression[mirna] < Q3:
            mediumExp.append(mirna)
        elif miRNAExpression[mirna] >= Q3:
            highExp.append(mirna)
    # verify that mirnas belong to all levels of expression
    assert len(lowExp) != 0
    assert len(moderateExp) != 0
    assert len(mediumExp) != 0
    assert len(highExp) != 0
    
    return [lowExp, moderateExp, mediumExp, highExp]



# use this function to create dict from miranda outputs
def SelectmiRNAsMirandaOutput(targetscan_seq_input_file, predicted_targets, expression_group, conservation, *seeds_mature):
    '''
    (file, file, list, str, *set, *file) -> dict
    Take the targetscan input sequence file, the miranda output, a list
    of mirnas sorted according to their expression level (low, moderate, medium or high),
    and return a dictionary with gene as key and a list with the number of
    predicted target sites, sequence length and number of target sites
    normalized by sequence length for all miRNAs or for conserved miRNAs only.
    If conservation = conserved, the a set of conserved seeds and the fasta file with 
    mature miRNAs is obtained from the optional arguments
    Precondition: files are generated to include valid or all chromos    
    '''
    
    
    # get seed set from optional parameter tuple
    if conservation == 'conserved':
        seeds = seeds_mature[0]
        # make a dict of mirna {name : seed} pairs
        mature_fasta = seeds_mature[1]
        mir_names = mirna_to_seed(mature_fasta)    
    
    # get the length of the sequences used to predict target sites {gene : seq_length}
    genes_length = get_domain_length_from_targetscan_input(targetscan_seq_input_file)
        
    # create a dict to store the target sites {gene: set(target1, target2, target3)} 
    target_counts = {}
    # create a dict with {gene : sequence length} pairs
    target_length = {}    
    
    # open file for reading
    infile = open(predicted_targets, 'r')
    # go through file
    for line in infile:
        if line.startswith('>>'):
            line = line.rstrip().split('\t')
            # get mirna, get rid of '>>' sign
            mirna = line[0][2:]
            # check if mirna in expression group
            if mirna in expression_group:
                # get rid of species code in mirna name
                mirna = mirna[4:]
                # get target gene
                gene = line[1]
                # get sequence length
                seq_length = int(line[8])
                # populate dict with gene : sequence length pairs
                if gene not in target_length:
                    target_length[gene] = seq_length
                # get positions
                positions = line[9].split()
                # check if consider all or conserved miRNA families
                if conservation == 'conserved':
                    # check that mirna seed in set of seeds
                    if mir_names[mirna] in seeds:
                        # record sites only if mirna seed in set of seeds
                        # populate dict
                        if gene not in target_counts:
                            # inititalise value with empty list
                            target_counts[gene] = set()
                        # add all sites to set
                        for pos in positions:
                            target_counts[gene].add(pos)
                elif conservation == 'all':
                    # record all sites
                    # populate dict
                    if gene not in target_counts:
                        # initialize value with empty set
                        target_counts[gene] = set()
                    for pos in positions:
                        target_counts[gene].add(pos)
                    
    # close file after reading
    infile.close()
    
    # create a dict to store the number of targets {gene: [N_sites, seq_length, N_sites/seq_length]}
    targets = {}
    
    # loop over genes with predicted targets in targetscan
    for gene in target_counts:
        # add number of target sites
        targets[gene] = [len(target_counts[gene])]
        # get the length of the sequence used to predict target sites
        seq_length = target_length[gene]
        # add sequence length to list
        targets[gene].append(seq_length)
        # add number of sites normalized by sequence length
        targets[gene].append(len(target_counts[gene]) / seq_length)
    
    # count 0 for genes that do not have any targets but that have a region > mininum length
    # loop over genes in targetscan seq input
    for gene in genes_length:
        # check that sites are not already recorded
        if gene not in targets:
            # targets were not predicted
            # check that sequence is greater than 6 bp
            if genes_length[gene] >= 7:
                # populate dict
                targets[gene] = [0, genes_length[gene], 0]

    return targets



# use this function to generate a score for target
def TargetScore(miRBaseFastaFile, NBins, species, miRBaseFile = 'miRNA.dat', ExpressionFile = 'mirna_read_count.txt'):
    '''
    (dict, dict, int, str) -> dict
    Take the dictionary with accession: mature names pairs for each species,
    the dictionary with mirna accession: expression pairs, the end boundary for 
    expression range and the species name and return a dict with mature name:
    score pairs. The score is simply the bin position of the mirna in the
    distribution of mirna expression.    
    '''
    
    # species is in the form Genus_species    
    
    # create a dictionary with {species: {mirna accession : mature names pairs}}
    AccessionNames = MatchmiRNAAccessionNumbers('accession', miRBaseFile)
    
    # create a dict {mature accession: mirna accession}
    miRNAAccessions = {}
    for accession in AccessionNames[species]:
        for mature in AccessionNames[species][accession]:
            miRNAAccessions[mature] = accession
    
    # match mature accessions to mature names {mature accession: mature names}
    MatureAccessions = MatchMatureAccessionsNames(miRBaseFastaFile)

    # create a dictionary {mirna accession: expression level}
    miRNAExpression = miRBAsemiRNAExpression(ExpressionFile)
    
    # create a dictionary with mature expression pairs {mature names: expression}
    MatureExpression = {}
    # loop over mature accessions
    for mature in MatureAccessions:
        # get mature name
        name = MatureAccessions[mature]
        # verify that mature has a mirna accession 
        assert mature in miRNAAccessions, 'mature accession should match a miRNA accession'
        # get mirna accession
        mirna = miRNAAccessions[mature]
        # get expression level
        if mirna in miRNAExpression:
            expression_level = miRNAExpression[mirna]
            # populate dict with mature name and expression level
            assert name not in MatureExpression, 'mature name should not already be recorded'
            MatureExpression[name] = expression_level
            
    # make a list of expression level 
    ExpressionLevel = [MatureExpression[mature] for mature in MatureExpression]
    ExpressionLevel.sort()
    
    # create a histogram with expression level 
    HistoCounts, HistoLimits = np.histogram(ExpressionLevel, range(0, NBins, 10))

    # Create a score based on miRNA expression to weight the importance of mirna targets
    # score is simply the bin position of the mirna expression
    Score = {}
    for mirna in MatureExpression:
        for i in range(0, len(HistoLimits) - 1):
            if MatureExpression[mirna] >= HistoLimits[i] and MatureExpression[mirna] < HistoLimits[i+1]:
                Score[mirna] = i+1
        if MatureExpression[mirna] >= HistoLimits[i+1]:
            Score[mirna] = len(HistoLimits)
    assert len(Score) == len(MatureExpression), 'scores are not recorded for some mirnas'

    return Score


# use this function to count the number of miranda target sites weighted by the mirna expression level
def WeightTargetsMirandaOutput(targetscan_seq_input_file, predicted_targets, Score):
    '''
    (file, file, dict) -> dict
    Take the targetscan input sequence file, the miranda output, the dictionary with 
    mirna: expression score pairs and return a dictionary with gene as key and
    a list with the number of predicted target sites, sequence length and number of target sites
    normalized by sequence length for all miRNAs or for conserved miRNAs only.
    All target counts are weighted by a score to take into account the expression of the cognate miRNA
    '''
    
    # get the length of the sequences used to predict target sites {gene : seq_length}
    genes_length = get_domain_length_from_targetscan_input(targetscan_seq_input_file)
        
    # create a dict to store the target sites {gene: number of weighted targets} 
    target_counts = {}
    # create a dict with {gene : sequence length} pairs
    target_length = {}    
    
    # open file for reading
    infile = open(predicted_targets, 'r')
    # go through file
    for line in infile:
        if line.startswith('>>'):
            line = line.rstrip().split('\t')
            # get mirna, get rid of '>>' sign
            mirna = line[0][2:]
            # get target gene
            gene = line[1]
            # get sequence length
            seq_length = int(line[8])
            # populate dict with gene : sequence length pairs
            if gene not in target_length:
                target_length[gene] = seq_length
            # get positions
            positions = line[9].split()
            # populate dict
            if gene not in target_counts:
                # initialize value 
                target_counts[gene] = 0
            # count the number of target sites weighted by the score of the mirna
            target_counts[gene] += (len(positions) * Score[mirna])
    # close file after reading
    infile.close()
    
    # create a dict to store the number of targets {gene: [N_sites, seq_length, N_sites/seq_length]}
    targets = {}
    # loop over genes with predicted targets
    for gene in target_counts:
        # add number of target sites
        targets[gene] = [target_counts[gene]]
        # get the length of the sequence used to predict target sites
        seq_length = target_length[gene]
        # add sequence length to list
        targets[gene].append(seq_length)
        # add number of sites normalized by sequence length
        targets[gene].append(target_counts[gene] / seq_length)
    
    # count 0 for genes that do not have any targets but that have a region > mininum length
    # loop over genes in targetscan seq input
    for gene in genes_length:
        # check that sites are not already recorded
        if gene not in targets:
            # targets were not predicted
            # check that sequence is greater than 6 bp
            if genes_length[gene] >= 7:
                # populate dict
                targets[gene] = [0, genes_length[gene], 0]

    return targets



# use this function to match mature accession to mature names for mirnas used to opredict sites
def MatchMatureAccessionsNames(miRBaseFastaFile):
    '''
    (file) -> dict
    Take the miRBase fasta file of unaligned mature sequences and return a dict
    matching the mature name with the mature accession
    '''
    
    # create a dict {accession: name}
    mature = {}
    # open file for reading
    infile = open(miRBaseFastaFile)
    for line in infile:
        if line.startswith('>'):
            line = line.rstrip().split()
            # get accession name
            accession = line[1]
            # get mirna name, remove '>'
            name = line[0][1:]
            mature[accession] = name
    infile.close()
    return mature
    
    