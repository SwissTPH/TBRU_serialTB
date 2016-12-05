#!/usr/bin/env python

#Clean up and re-format the data into DataFrame

#NB written for Python 3.5+

################
# DEPENDENCIES #
################

import os
import pickle
import pandas as pd
import numpy as np
import serial_functions as serf #custom functions (./serial_functions.py)

#############################
#  LOAD PRE-PROCESSED DATA  #
#############################

CWD = os.getcwd() #get current working directory

INTERIM_PATH = os.path.abspath(os.path.join(CWD, os.pardir, os.pardir, 'data', 'interim')) #path to interim data

#load individual data
PATIENT_DATA = pickle.load(open('{}{}2_patient_data.pkl'.format(INTERIM_PATH,os.sep),'rb'))
ANNOTATION_DICT = pickle.load(open('{}{}2_unfixed_annotation.pkl'.format(INTERIM_PATH,os.sep),'rb'))
DATA_DICT = pickle.load(open('{}{}2_unfixed_mutations.pkl'.format(INTERIM_PATH,os.sep),'rb'))
FIXED_ANNOTATION_DICT = pickle.load(open('{}{}2_fixed_annotation.pkl'.format(INTERIM_PATH,os.sep),'rb'))
FIXED_ARRAY = pickle.load(open('{}{}2_fixed_mutations.pkl'.format(INTERIM_PATH,os.sep),'rb'))

###############################
#  CHRUN PRE-PROCESSED vSNPs  #
###############################

TEMP_LIST = []

_TN = {'nonsynonymous': 'NSY', 'synonymous': 'SYN'}

for i,v in enumerate(PATIENT_DATA['PATIENTS']):
    for spot in range(1,8):
        #patient treated efficaciously
        _tp = int(v in ['Patient09', 'Patient10', 'Patient11', 'Patient12'])
        #get all alleles that have signal at a given timepoint
        _data = PATIENT_DATA['UNFIXED_ARRAYS'][i][:,spot][PATIENT_DATA['UNFIXED_ARRAYS'][i][:,spot]>0]
        _loci = PATIENT_DATA['UNFIXED_ARRAYS'][i][:,0][PATIENT_DATA['UNFIXED_ARRAYS'][i][:,spot]>0]
        if len(_data)>0:
            for _frequency,_locus in zip(_data,_loci):
                _locus = str(int(_locus))
                #annotate the SNP
                _annotation = ANNOTATION_DICT[_locus][0]
                #Determine the wild type codon if applicable, otherwise return -
                _wt_codon = '-'
                if _annotation[1] in ['synonymous','nonsynonymous']:
                    _wt_codon = _annotation[4].split('-')[0]
                #Collate all the data into a temporary list
                TEMP_LIST.append([v, PATIENT_DATA['TM'][spot-1], PATIENT_DATA['TP2'][i],
                             _frequency, np.log(_frequency), _tp,
                             _locus, _annotation[0],
                             _TN.get(_annotation[1], _annotation[1]),
                             _annotation[-1], _wt_codon
                             ])

##########################
#  CONVERT TO DATAFRAME  #
##########################

ALL_TEMP = pd.DataFrame({'PATIENT_ID': [x[0] for x in TEMP_LIST],
                    'TIME': [x[1] for x in TEMP_LIST],
                    'RESISTANCE': [x[2] for x in TEMP_LIST],
                    'FREQUENCY': [x[3] for x in TEMP_LIST],
                    'LOG_FREQ': [x[4] for x in TEMP_LIST],
                    'NON_EFFICACIOUS': [x[5] for x in TEMP_LIST],
                    'LOCUS': [x[6] for x in TEMP_LIST],
                    'GENE': [x[7] for x in TEMP_LIST],
                    'SNP_TYPE': [x[8] for x in TEMP_LIST],
                    'TRANSITION': [x[9] for x in TEMP_LIST],
                    'WT_CODON': [x[10] for x in TEMP_LIST]
                    })

#Tag alleles
ALL_TEMP['ALLELE_TAG'] = ['{}_{}'.format(x,y) for (x,y) in zip(ALL_TEMP.PATIENT_ID, ALL_TEMP.LOCUS)]

#Identify vSNPs with a frequency above 1.5%
LOCI_ABOVE_CUTOFF = list(set(ALL_TEMP.ALLELE_TAG[ALL_TEMP.FREQUENCY>0.015]))

#Identify loci that have a frequency above 1.5% in at least one isolate from a patient
ALL_TEMP['KEEP'] = [int(x in LOCI_ABOVE_CUTOFF) for x in ALL_TEMP.ALLELE_TAG]

#Filter away all vSNPs that never appear above the 1.5% frequency threshold
ALL = ALL_TEMP[(ALL_TEMP.KEEP==1)].copy()
ALL.drop('KEEP', axis=1, inplace=True) #delete the 'KEEP' column
ALL.drop('ALLELE_TAG', axis=1, inplace=True) #delete the 'ALLELE_TAG' column
ALL.reset_index(drop=True, inplace=True) #reset the index counter

ALL.to_csv('{}3_ALLELE_data.csv'.format(INTERIM_PATH))

###############################
#  CHRUN PRE-PROCESSED fSNPs  #
###############################

TEMP_LIST = []

for i,v in enumerate(PATIENT_DATA['PATIENTS']):
    #patient treated efficaciously
    _tp = int(v in ['Patient09', 'Patient10', 'Patient11', 'Patient12'])
    #get all fSNPs for a patient
    _loci = np.unique(FIXED_ARRAY[FIXED_ARRAY[:,i+1]>0][:,0])
    for _locus in _loci:
        _locus = str(int(_locus))
        #annotate the SNP
        _annotation = FIXED_ANNOTATION_DICT[_locus][0]
        #Determine the wild type codon if applicable, otherwise return '-'
        _wt_codon = '-'
        if _annotation[1] in ['synonymous','nonsynonymous']:
            _wt_codon = _annotation[4].split('-')[0]
        #Collate all the data into a temporary list
        TEMP_LIST.append([v, PATIENT_DATA['TP2'][i],
                          _tp, _locus, _annotation[0],
                          _TN.get(_annotation[1], _annotation[1]),
                          _annotation[-1], _wt_codon
                          ])

##########################
#  CONVERT TO DATAFRAME  #
##########################

FIX = pd.DataFrame({'PATIENT_ID': [x[0] for x in TEMP_LIST],
                    'RESISTANCE': [x[1] for x in TEMP_LIST],
                    'NON_EFFICACIOUS': [x[2] for x in TEMP_LIST],
                    'LOCUS': [x[3] for x in TEMP_LIST],
                    'GENE': [x[4] for x in TEMP_LIST],
                    'SNP_TYPE': [x[5] for x in TEMP_LIST],
                    'TRANSITION': [x[6] for x in TEMP_LIST],
                    'WT_CODON': [x[7] for x in TEMP_LIST]
                    })

FIX.to_csv('{}3_FIXED_data.csv'.format(INTERIM_PATH))
