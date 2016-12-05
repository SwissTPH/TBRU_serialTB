#!/usr/bin/env python

#Parse ART simulation data

#ART version used:

#ART command used to for running simulations
#

#500 iterations were run and passed through the pipeline

##############
#KERNEL SETUP#
##############

import pickle
import pandas as pd
import datetime

###############
#PATHS TO DATA#
###############

PATH_TO_RESULT = '../../data/raw/ART/'
ART_LIST = '../../data/interim/ART_files.list'

################
#OPEN AND PARSE#
################

INDEX_KEY = 0
ART_DATA = {} #data will be stored as a dictionary: easy conversion to DataFrame

for result in open(ART_LIST):
    #PARSE FINAL RESULT "X.mixupf.lablef"
    result = result.strip()
    _SIM_ID = result.split('.')[0]
    for line in open('{}{}'.format(PATH_TO_RESULT,result)):
        line = line.strip()
        split = line.split()

        _locus = int(split[3])
        _ref_base = line.split()[4][1:-1]
        _alt_base = line.split()[5][1:-1]
        _freq = float(line.split()[6][1:-2])/100
        _coverage = int(split[7])

        ART_DATA[INDEX_KEY] = {'SIMULATION': _SIM_ID, 'LOCUS': _locus,
                               'REF_BASE': _ref_base, 'ALT_BASE': _alt_base,
                               'FREQUENCY': round(_freq,4), 'STAGE': 'FINAL',
                               'COVERAGE': _coverage}
        INDEX_KEY+=1

    #PARSE FILTERED RESULT "X.mixupf"
    result = result[:-7]
    for line in open('{}{}'.format(PATH_TO_RESULT,result)):
        line = line.strip()
        split = line.split()

        _locus = int(split[1])
        _ref_base = line.split()[2]
        _alt_base = line.split()[3]
        _freq = float(line.split()[4][:-1])/100
        _coverage = int(split[5])

        ART_DATA[INDEX_KEY] = {'SIMULATION': _SIM_ID, 'LOCUS': _locus,
                               'REF_BASE': _ref_base, 'ALT_BASE': _alt_base,
                               'FREQUENCY': round(_freq,4), 'STAGE': 'FILTER',
                               'COVERAGE': _coverage}
        INDEX_KEY+=1

    #PARSE BASIC OUTPUT OF VARIANT CALLING "X.mixup"
    result = result[:-1]
    for line in open('{}{}'.format(PATH_TO_RESULT,result)):
        line = line.strip()
        split = line.split()

        _locus = int(split[1])
        _ref_base = line.split()[2]
        _alt_base = line.split()[3]
        _freq = float(line.split()[4][:-1])/100
        _coverage = int(split[5])

        ART_DATA[INDEX_KEY] = {'SIMULATION': _SIM_ID, 'LOCUS': _locus,
                               'REF_BASE': _ref_base, 'ALT_BASE': _alt_base,
                               'FREQUENCY': round(_freq,4), 'STAGE': 'BASIC',
                               'COVERAGE': _coverage}
        INDEX_KEY+=1



ART_DF = pd.DataFrame(ART_DATA).T

NOW = datetime.datetime.now()
#TARGET = '../../data/interim/{}_ART_simulation'.format(NOW.strftime('%y%m%d'))
TARGET = '../../data/interim/1_ART_simulation'

pickle.dump(ART_DATA, open('{}.pkl'.format(TARGET),'wb'))
ART_DF.to_csv('{}.csv'.format(TARGET))
