#!/usr/bin/env python

#Parse FINAL, annotated output of patient vSNP data

#These vSNP calls emerge from the first part of the pipeline
#they are filtered by:
#        - base-calling quality,
#        - mapping quality,
#        - minimal overall coverage,
#        - at least 2 forward and 2 reverse reads,
#        - strand filter,
#        - minimum frequency: 0.5%
#        - remove problematic loci,
#        - exclude calls with read-end bias,
#        - annotated to H37Rv, using H37Rv loci
#          if a vSNP locus was absent in H37Rv it was given a % and reported
#          as it was in CCDC5079.

#12 patients with several samples each

##############
#KERNEL SETUP#
##############

import os
import pickle
import pandas as pd
import datetime

###############
#PATHS TO DATA#
###############

CWD = os.getcwd() #get current working directory


INTERIM_PATH = os.path.abspath(os.path.join(CWD, os.pardir, os.pardir, 'data', 'interim'))
PATH_TO_RESULT = os.path.abspath(os.path.join(CWD, os.pardir, os.pardir, 'data', 'raw', 'FINAL_vSNP'))

FINAL_LIST = '{}{}FINAL.list'.format(INTERIM_PATH,os.sep)
SNPS = '{}{}SNP_COVERAGE.pkl'.format(INTERIM_PATH,os.sep)

################
#OPEN AND PARSE#
################

INDEX_KEY = 0
FINAL_DATA = {} #data will be stored as a dictionary: easy conversion to DataFrame


try:
    SNP_COVERAGE = pickle.load(open(SNPS,'rb'))
except IOError:
    print('Cannot open - make sure the following file exists:', SNPS)

for result in open(FINAL_LIST):
    #PARSE BASIC RESULT "X.snp"
    result = result.strip()
    _SAMPLE_ID = result.split('.')[0]
    _PATIENT_ID = _SAMPLE_ID[:2]
    _TIME_ID = _SAMPLE_ID[2:]
    for line in open('{}{}'.format(PATH_TO_RESULT,result)):
        line = line.strip()

        split = line.split()

        _locus = int(split[0])
        _ref_base = line.split()[1]
        _alt_base = line.split()[2]
        _freq = float(split[3][:-1])/100
        #Get the coverage for a SNP from a specific sample (nested dictionary)
        try:
            _coverage = SNP_COVERAGE.get(_locus,None).get(_SAMPLE_ID, None)
        except:
            _coverage = None


        FINAL_DATA[INDEX_KEY] = {'PATIENT': _PATIENT_ID, 'LOCUS': _locus,
                                 'REF_BASE': _ref_base, 'ALT_BASE': _alt_base,
                                 'FREQUENCY': round(_freq,4), 'STAGE': 'FINAL',
                                 'COVERAGE': _coverage, 'TIME': _TIME_ID,
                                 'SAMPLE': _SAMPLE_ID}
        INDEX_KEY+=1


FINAL_DF = pd.DataFrame(FINAL_DATA).T

NOW = datetime.datetime.now()
#TARGET = '../../data/interim/{}_FINAL'.format(NOW.strftime('%y%m%d'))
TARGET = '{}{}1_FINAL'.format(INTERIM_PATH,os.sep)

pickle.dump(FINAL_DATA, open('{}.pkl'.format(TARGET),'wb'))
FINAL_DF.to_csv('{}.csv'.format(TARGET))
