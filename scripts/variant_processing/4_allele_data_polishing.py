#!/usr/bin/env python

#Final processing of the Allele data prior to analysis

#NB written for Python 3.5+

################
# DEPENDENCIES #
################

import os
import pickle
import pandas as pd
import numpy as np

#############################
#  LOAD PRE-PROCESSED DATA  #
#############################

CWD = os.getcwd() #get current working directory

PROCESSED_PATH = os.path.abspath(os.path.join(CWD, os.pardir, os.pardir, 'data', 'processed'))
INTERIM_PATH = os.path.abspath(os.path.join(CWD, os.pardir, os.pardir, 'data', 'interim'))

ALL_TEMP = pd.read_csv('{}{}3_ALLELE_data.csv'.format(INTERIM_PATH, os.sep), index_col=0)

#Remove manually identified noisy vSNP:
# 2347614, 2347615, 2637298, 2637459
LOCI_TO_EXCLUDE = [2347614, 2347615, 2637298, 2637459]

ALL_TEMP['EXCLUDE'] = [int(x in LOCI_TO_EXCLUDE) for x in ALL_TEMP.LOCUS]

#Because only a handful of patients have data for 8 weeks+, filter data by TIME
#Remove vSNPs at LOCI_TO_EXCLUDE
ALL = ALL_TEMP[(ALL_TEMP.EXCLUDE==0)&(ALL_TEMP.TIME<12)].copy()
ALL.drop('EXCLUDE', axis=1, inplace=True) #delete the 'EXCLUDE' column

#Save the resulting dataframe
ALL.to_csv('{}{}ALLELE_data.csv'.format(PROCESSED_PATH, os.sep))


##########################
#  UPDATE ND.ARRAY DATA  #
##########################

PATIENT_DICT = pickle.load(open('{}{}2_patient_data.pkl'.format(INTERIM_PATH, os.sep),'rb'))

#Make the nd.array with vSNPs consistent with the DataFrame
for ind,patient in enumerate(PATIENT_DICT['PATIENTS']):
    #get all the loci in float format
    _alleles_to_keep = np.array(ALL.LOCUS[ALL.PATIENT_ID==patient], dtype=float)
    _alleles_to_filter = PATIENT_DICT['UNFIXED_ARRAYS'][ind][:,0]
    #determine which data to keep
    _filter_indices = \
    np.array([x in _alleles_to_keep for x in _alleles_to_filter])
    #filter the arrays in the dictionary including the timepoint data
    PATIENT_DICT['UNFIXED_ARRAYS'][ind] = \
    _alleles_to_filter = PATIENT_DICT['UNFIXED_ARRAYS'][ind][_filter_indices,:-2]

pickle.dump(PATIENT_DICT, open('{}{}ALLELE_info.pkl'.format(PROCESSED_PATH, os.sep),'wb'))
