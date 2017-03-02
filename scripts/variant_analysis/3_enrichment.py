#!/usr/bin/env python

#Estimate whether excessive mutations have been

#NB written for Python 3.5+

################
# DEPENDENCIES #
################

import os
import pickle
import pandas as pd
import numpy as np
import scipy.stats as ss
import serial_functions as serf #custom functions (./serial_functions.py)

#############################
#  LOAD PRE-PROCESSED DATA  #
#############################

PROCESSED_PATH = os.path.abspath(os.path.join(CWD, os.pardir, os.pardir, 'data', 'processed'))
INTERIM_PATH = os.path.abspath(os.path.join(CWD, os.pardir, os.pardir, 'data', 'interim'))

ALL = pd.read_csv('{}ALLELE_data.csv'.format(PROCESSED_PATH), index_col=0)
EGS = pd.read_csv('{}PREDEFINED_gene_sets.csv'.format(EXTERNAL_PATH), index_col=0)

#Prepare the necessary resources from the MTBC genome
REF_ANNOTATIONS = {}
for line in open('{}H37Rv_annotation2sytems.ptt'.format(EXTERNAL_PATH)):
    split = line.strip().split('\t')
    if len(split)<=6 and len(split)>3 and split[0][0].isdigit():
        REF_ANNOTATIONS[int(split[1])-1] = [int(split[2]), split[3], split[4]]
    if len(split)>6 and split[0][0].isdigit():
        REF_ANNOTATIONS[int(split[1])-1] = [int(split[2]), split[3], split[7]]

LOCUS_MATRIX = np.array([[k,v[0]] for (k,v) in REF_ANNOTATIONS.items()], dtype=float)
LOCUS_MATRIX = LOCUS_MATRIX[np.argsort(LOCUS_MATRIX[:,0])]

LOCI = np.array([REF_ANNOTATIONS[int(x)][-1] for x in LOCUS_MATRIX[:,0]])
SIZE_DICT = {x:float(LOCUS_MATRIX[ind][1])-LOCUS_MATRIX[ind][0] for (ind,x) in enumerate(LOCI)}
MTB_GENOME = ''.join([line.strip() for line in open('{}MTB_anc.fa'.format(EXTERNAL_PATH)) if line[0]!='>'])

#Define the excluded geneset (Coscolla et al, 2015)
excluded = list(EGS.GENE[EGS.GENE_SET=='Excluded'])

#Define the drug resistance associated gene set (Farhat et al, 2013 "approach overlap genes"
# and Zhang et al 2013 - "Drug resistance group 1")
DR_set = list(EGS.GENE[(EGS.GENE_SET=='DR_associated')&(EGS.EXCLUDED==0)])

#Define the T-cell antigen gene set (Coscolla et al, 2015)
AG_set = list(EGS.GENE[(EGS.GENE_SET=='Antigens')&(EGS.EXCLUDED==0)])

#Define the drug resistance gene set (Walker et al, 2015)
Walker_DR = list(EGS.GENE[(EGS.GENE_SET=='DR_genes')&(EGS.EXCLUDED==0)])

#Define Mycolate superpathway genes (O'Neill et al, 2015)
ONeillMycolate = list(EGS.GENE[(EGS.GENE_SET=='Mycolate_superpathway')&(EGS.EXCLUDED==0)])

###########################
#  DRUG RESISTANCE GENES  #
###########################

DR_output_E = serf.excess_mutation(Walker_DR, ALL[(ALL.NON_EFFICACIOUS==0)],
                               len(MTB_GENOME), SIZE_DICT, exclude_check=excluded)
DR_output_NE = serf.excess_mutation(Walker_DR, ALL[(ALL.NON_EFFICACIOUS==1)],
                               len(MTB_GENOME), SIZE_DICT, exclude_check=excluded)

print('WALKER DR GENESET\n')
print('OOOO\tEfficaciously treated\tOOOO')
print('--- COUNTS ---')
print('SNPS AT DIAGNOSIS:\n%s\n' % DR_output_E['Data'][0])
print('SNPS DURING TREATMENT:\n%s\n' % DR_output_E['Data'][1])
print('--- STATISTICS ---')
print('Excess binomial: %.3f' % DR_output_E['Excess_binomial'])
print('Excess NSY: %.3f\n' % DR_output_E['NS_binomial'])

print('OOOO\tNon-efficaciously treated\tOOOO')
print('--- COUNTS ---')
print('SNPS AT DIAGNOSIS:\n%s\n' % DR_output_NE['Data'][0])
print('SNPS DURING TREATMENT:\n%s\n' % DR_output_NE['Data'][1])
print('--- STATISTICS ---')
print('Excess binomial: %.3f' % DR_output_NE['Excess_binomial'])
print('Excess NSY: %.3f\n' % DR_output_NE['NS_binomial'])

######################################
#  DRUG RESISTANCE-ASSOCIATED GENES  #
######################################

DR_output_E = serf.excess_mutation(DR_set, ALL[(ALL.NON_EFFICACIOUS==0)],
                               len(MTB_GENOME), SIZE_DICT, exclude_check=excluded)
DR_output_NE = serf.excess_mutation(DR_set, ALL[(ALL.NON_EFFICACIOUS==1)],
                               len(MTB_GENOME), SIZE_DICT, exclude_check=excluded)

print('DR-associated GENESET\n')
print('OOOO\tEfficaciously treated\tOOOO')
print('--- COUNTS ---')
print('SNPS AT DIAGNOSIS:\n%s\n' % DR_output_E['Data'][0])
print('SNPS DURING TREATMENT:\n%s\n' % DR_output_E['Data'][1])
print('--- STATISTICS ---')
print('Excess binomial: %.3f' % DR_output_E['Excess_binomial'])
print('Excess NSY: %.3f\n' % DR_output_E['NS_binomial'])

print('OOOO\tNon-efficaciously treated\tOOOO')
print('--- COUNTS ---')
print('SNPS AT DIAGNOSIS:\n%s\n' % DR_output_NE['Data'][0])
print('SNPS DURING TREATMENT:\n%s\n' % DR_output_NE['Data'][1])
print('--- STATISTICS ---')
print('Excess binomial: %.3f' % DR_output_NE['Excess_binomial'])
print('Excess NSY: %.3f\n' % DR_output_NE['NS_binomial'])

###########################
#  Mycolate superpathway  #
###########################

output_E = serf.excess_mutation(ONeillMycolate, ALL[(ALL.NON_EFFICACIOUS==0)],
                               len(MTB_GENOME), SIZE_DICT, exclude_check=excluded)
output_NE = serf.excess_mutation(ONeillMycolate, ALL[(ALL.NON_EFFICACIOUS==1)],
                               len(MTB_GENOME), SIZE_DICT, exclude_check=excluded)

print('Mycolate superpathway\n')
print('OOOO\tEfficaciously treated\tOOOO')
print('--- COUNTS ---')
print('SNPS AT DIAGNOSIS:\n%s\n' % output_E['Data'][0])
print('SNPS DURING TREATMENT:\n%s\n' % output_E['Data'][1])
print('--- STATISTICS ---')
print('Excess binomial: %.3f' % output_E['Excess_binomial'])
print('Excess NSY: %.3f\n' % output_E['NS_binomial'])

print('OOOO\tNon-efficaciously treated\tOOOO')
print('--- COUNTS ---')
print('SNPS AT DIAGNOSIS:\n%s\n' % output_NE['Data'][0])
print('SNPS DURING TREATMENT:\n%s\n' % output_NE['Data'][1])
print('--- STATISTICS ---')
print('Excess binomial: %.3f' % output_NE['Excess_binomial'])
print('Excess NSY: %.3f\n' % output_NE['NS_binomial'])

#####################
#  T-cell antigens  #
#####################

output_E = serf.excess_mutation(AG_set, ALL[(ALL.NON_EFFICACIOUS==0)],
                               len(MTB_GENOME), SIZE_DICT, exclude_check=excluded)
output_NE = serf.excess_mutation(AG_set, ALL[(ALL.NON_EFFICACIOUS==1)],
                               len(MTB_GENOME), SIZE_DICT, exclude_check=excluded)

print('T-cell antigens\n')
print('OOOO\tEfficaciously treated\tOOOO')
print('--- COUNTS ---')
print('SNPS AT DIAGNOSIS:\n%s\n' % output_E['Data'][0])
print('SNPS DURING TREATMENT:\n%s\n' % output_E['Data'][1])
print('--- STATISTICS ---')
print('Excess binomial: %.3f' % output_E['Excess_binomial'])
print('Excess NSY: %.3f\n' % output_E['NS_binomial'])

print('OOOO\tNon-efficaciously treated\tOOOO')
print('--- COUNTS ---')
print('SNPS AT DIAGNOSIS:\n%s\n' % output_NE['Data'][0])
print('SNPS DURING TREATMENT:\n%s\n' % output_NE['Data'][1])
print('--- STATISTICS ---')
print('Excess binomial: %.3f' % output_NE['Excess_binomial'])
print('Excess NSY: %.3f\n' % output_NE['NS_binomial'])
