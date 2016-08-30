# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 12:15:31 2016

@author: RJovelin
"""



# use this script to download chromosome sequence files and unzip them

import os

# counter number of files
total = 0
infile = open('GRCH37_FTP_links.txt')
# loop over file, get the ftp link, download the corresponding file
for line in infile:
    line = line.rstrip()
    if line != '':
        total += 1
        os.system('wget ' + line)
infile.close()

# make a list of gz files
files = [i for i in os.listdir() if i[-3:] == '.gz']
print(len(files))
assert len(files) == total, 'numbers of files to unzip and to download do not match'
