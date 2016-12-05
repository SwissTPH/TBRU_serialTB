#!/usr/bin/env python

#Parse data from the sequencing of individual MTB colonies

#Two MTB colonies were grown in 7H9 and sequenced to great depth twice.
#This yielded four sequencing samples that we plan to use for determining
#the need for an empirical variant frequency cutoff.

#Meta data on sequencing...

##############
#KERNEL SETUP#
##############

import pickle
import pandas as pd
import datetime

###############
#PATHS TO DATA#
###############

PATH_TO_RESULT = '../../data/raw/single_colony/'
COLONY_LIST = '../../data/interim/single_colony_files.list'

################
#OPEN AND PARSE#
################

INDEX_KEY = 0
COLONY_DATA = {} #data will be stored as a dictionary: easy conversion to DataFrame

for result in open(COLONY_LIST):
    #PARSE FINAL OUTPUT OF PIPELINE "x.mixupf.lable"
    result = result.strip()
    _COLONY_ID = result.split('.')[0]
    for line in open('{}{}'.format(PATH_TO_RESULT,result)):
        line = line.strip()
        split = line.split()

        _locus = int(split[3])
        _ref_base = line.split()[4][1:-1]
        _alt_base = line.split()[5][1:-1]
        _freq = float(line.split()[6][1:-2])/100
        _coverage = int(split[7])

        COLONY_DATA[INDEX_KEY] = {'COLONY': _COLONY_ID, 'LOCUS': _locus,
                                  'REF_BASE': _ref_base, 'ALT_BASE': _alt_base,
                                  'FREQUENCY': round(_freq,4), 'STAGE': 'FINAL',
                                  'COVERAGE': _coverage}
        INDEX_KEY+=1

    #PARSE BASIC OUTPUT OF VARIANT CALLING "x.mixup"
    result = result[:-7]
    for line in open('{}{}'.format(PATH_TO_RESULT,result)):
        line = line.strip()
        split = line.split()

        _locus = int(split[1])
        _ref_base = line.split()[2]
        _alt_base = line.split()[3]
        _freq = float(line.split()[4][:-1])/100
        _coverage = int(split[5])


        COLONY_DATA[INDEX_KEY] = {'COLONY': _COLONY_ID, 'LOCUS': _locus,
                                  'REF_BASE': _ref_base, 'ALT_BASE': _alt_base,
                                  'FREQUENCY': round(_freq,4), 'STAGE': 'BASIC',
                                  'COVERAGE': _coverage}
        INDEX_KEY+=1

COLONY_DF = pd.DataFrame(COLONY_DATA).T

NOW = datetime.datetime.now()
#TARGET = '../../data/interim/{}_colony_sequencing'.format(NOW.strftime('%y%m%d'))
TARGET = '../../data/interim/1_colony_sequencing'

pickle.dump(COLONY_DATA, open('{}.pkl'.format(TARGET),'wb'))
COLONY_DF.to_csv('{}.csv'.format(TARGET))
