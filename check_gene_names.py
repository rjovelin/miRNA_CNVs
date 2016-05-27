# -*- coding: utf-8 -*-
"""
Created on Sat Oct  3 17:10:27 2015

@author: Richard
"""

import os
from CNV_miRNAs import *

files = [i for i in os.listdir() if 'gff3' in i]


for filename in files:
    print(filename)
    longest = get_longest_mRNA(filename)
    print('longest ', len(longest))
    gene_names = geneID_to_gene_name(filename)
    print('gene id to name ', len(gene_names))
    rna_to_gene = match_mRNA_to_gene(filename)
    print('rnas to genes ', len(rna_to_gene))
    gene_ids = set(rna_to_gene[rna] for rna in longest)
    print('gene_ids ', len(gene_ids))
    print('dupli gene ids in longest ', len(longest) - len(gene_ids))
    names = set(gene_names[ID] for ID in gene_ids)
    print('names ', len(names))
    print('dupli names ', len(gene_ids) - len(names))