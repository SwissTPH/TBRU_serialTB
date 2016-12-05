#!/usr/bin/env python

#Script written by Andrej Trauner, 30/07/2015
# The script is to be used on $.merge.lable files to emerge from the SNP calling pipeline assembled by Qingyun.

# Usage: python 5_lable_filter.py $.merge.lable list.omit

#################
#    ARGUMENTS  #
#################

# $.merge.lable: SNP calls to be filtered
# list.omit: text file containing all the sites to be omitted.


###################
#    DESCRIPTION  #
###################

# 1. filter out tail SNPs
# 2. filter out SNPs with fewer than 0.5% frequencies
# 3. filter out all the PPE/INTEGRASE... based on the positions I got from Qingyun

##############
#    OUTPUT  #
##############

# $.filt pipeline SNPs filtered for abundance, call distribution and call location.

import sys
from collections import Counter


data_file = sys.argv[1]
prefix = data_file.split('.')[0] #get the file prefix

positions_to_exclude = sys.argv[2]

min_frequency = 0.5

#Generate the omission list

omit_these = [line.strip() for line in open(positions_to_exclude)] #strip command gets rid of empty space around the value

#Set up the output file:

output_filename = prefix+'.filter.snp'

output_file = open(output_filename,'w')

#Run through and filter:

for line in open(data_file):
    split = line.strip().split(' ')
    #SNP valid: split[1], SNP position[3], REFERENCE base: split[4]
    #Base calls: split[6], Percentage: split[9]
    freq = float(line.split(' ')[9][1:-2])
    alternative_calls = [
        (Counter(split[6].upper()).get(x, 0), x) #count all symbols in the string
        for x in ['A', 'C', 'G', 'T']] #circle through all bases
    alt_base = sorted(alternative_calls)[-1][1] #keep the most frequent alternative call.
    try:
        ref_base = split[4][1] #this is here in case the base is in qutations e.g. "C"
    except:
        ref_base = split[4]
    if split[1]=='1' and freq>min_frequency and split[3] not in omit_these:
        entry = '{}\t{}\t{}\n'.format(split[3], ref_base, alt_base)
        output_file.write(entry)

output_file.close()
