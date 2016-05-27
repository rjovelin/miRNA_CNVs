# -*- coding: utf-8 -*-
"""
Created on Sun Sep 27 00:35:34 2015

@author: Richard
"""


import os


# make a list of sequence files
seq_files = ['bta_seq_files.txt', 'cfa_seq_files.txt',
             'gga_seq_files.txt', 'hsa_seq_files.txt',
             'mml_seq_files.txt', 'mmu_seq_files.txt',
             'ptr_seq_files.txt', 'rno_seq_files.txt']

# loop over list of file
for filename in seq_files:
    # open file for reading
    infile = open(filename, 'r')
    # get the species name
    species = infile.readline().rstrip()
    # get the ftp directory
    ftp_dir = infile.readline().rstrip()
    # make a list of fasta files
    fasta_files = []
    for line in infile:
        line = line.rstrip()
        if line != '':
            fasta_files.append(line)
    # close file after reading
    infile.close()
    # create a directory with species name
    os.system('mkdir ' + species)
    # change directory
    os.chdir('./' + species + '/')
    # download all files
    for i in fasta_files:
        fetch_command = 'wget ' + ftp_dir + i
        print(fetch_command)
        os.system(fetch_command)
    # unzip files in current directory
    for zipfile in os.listdir():
        os.system('gunzip ' + zipfile)
    # go back to parent directory
    os.chdir('../')
        
# download GFF files
infile = open('GFF_files.txt', 'r')
# loop over file
for line in infile:
    line = line.rstrip()
    if line != '':
        # download file
        os.system('wget ' + line)
# close file
infile.close()

# unzip gff files
gff_files = [i for i in os.listdir() if 'gff3' in i]
for i in gff_files:
    os.system('gunzip ' + i)
    
print('downloaded all files')