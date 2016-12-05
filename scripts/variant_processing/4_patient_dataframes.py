#!/usr/bin/env python

#Processing Allele data per timepoint and patient

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

PROCESSED_PATH = os.path.abspath(os.path.join(CWD, os.pardir, os.pardir, 'data', 'processed'))
INTERIM_PATH = os.path.abspath(os.path.join(CWD, os.pardir, os.pardir, 'data', 'interim'))

ALL = pd.read_csv('{}{}ALLELE_data.csv'.format(PROCESSED_PATH, os.sep), index_col=0)
PATIENT_DATA = pickle.load(open('{}{}ALLELE_info.pkl'.format(PROCESSED_PATH, os.sep),'rb'))
ANNOTATION_DICT = pickle.load(open('{}{}2_unfixed_annotation.pkl'.format(INTERIM_PATH, os.sep),'rb'))

CCDC5079_genome = 4414325 #based on NC_021251 (CCDC5079)

#####################################
#  PERFORM TIME_POINT CALCULATIONS  #
#####################################

TIMEPOINT_DATA = {'PATIENT_ID': [],
                  'TIME': [],
                  'MEAN_FREQUENCY': [],
                  'LOG_MEAN_FREQUENCY': [],
                  'MEDIAN_FREQUENCY': [],
                  'VARIANCE_FREQUENCY': [],
                  'vSNP_COUNT': [],
                  'HETEROZYGOSITY': [],
                  'DEPTH': [],
                  'MEAN_SYN_FREQUENCY': [],
                  'MEAN_NSY_FREQUENCY': [],
                  'pNS_vSNP': [],
                  'NEUTRAL_pNS_vSNP': [],
                  'NON_EFFICACIOUS': []
                  }

for _patient in PATIENT_DATA['PATIENTS']:
    _timepoints = list(set(ALL.TIME[ALL.PATIENT_ID==_patient]))
    for _timepoint in _timepoints:
        _frequencies = list(ALL.FREQUENCY[(ALL.TIME==_timepoint)&(ALL.PATIENT_ID==_patient)])
        _mean_frequency = np.mean(_frequencies)
        _median_frequency = np.median(_frequencies)
        _variance = np.var(_frequencies)
        _heterozygosity = serf.heterozygosity(np.array(_frequencies),
                                              CCDC5079_genome)
        _effective = int(_patient in ['Patient09', 'Patient10', 'Patient11', 'Patient12'])
        _NSY_mean = np.mean(ALL.FREQUENCY[(ALL.TIME==_timepoint)&(ALL.PATIENT_ID==_patient)&(ALL.SNP_TYPE=='NSY')])
        _SYN_mean = np.mean(ALL.FREQUENCY[(ALL.TIME==_timepoint)&(ALL.PATIENT_ID==_patient)&(ALL.SNP_TYPE=='SYN')])

        TIMEPOINT_DATA['PATIENT_ID'].append(_patient)
        TIMEPOINT_DATA['TIME'].append(_timepoint)
        TIMEPOINT_DATA['MEAN_FREQUENCY'].append(_mean_frequency)
        TIMEPOINT_DATA['LOG_MEAN_FREQUENCY'].append(np.log(_mean_frequency))
        TIMEPOINT_DATA['MEDIAN_FREQUENCY'].append(_median_frequency)
        TIMEPOINT_DATA['VARIANCE_FREQUENCY'].append(_variance)
        TIMEPOINT_DATA['MEAN_SYN_FREQUENCY'].append(_SYN_mean)
        TIMEPOINT_DATA['MEAN_NSY_FREQUENCY'].append(_NSY_mean)
        TIMEPOINT_DATA['vSNP_COUNT'].append(len(_frequencies))
        TIMEPOINT_DATA['HETEROZYGOSITY'].append(_heterozygosity)
        TIMEPOINT_DATA['DEPTH'].append(PATIENT_DATA['DEPTH'][_patient][_timepoint])
        TIMEPOINT_DATA['NON_EFFICACIOUS'].append(_effective)
        TIMEPOINT_DATA['pNS_vSNP'].append(_effective)
        TIMEPOINT_DATA['NEUTRAL_pNS_vSNP'].append(_effective)

TIMEPOINT_DF = pd.DataFrame(TIMEPOINT_DATA)

TIMEPOINT_DF.to_csv('{}{}TIMEPOINT_data.csv'.format(PROCESSED_PATH, os.sep))

##################################
#  PERFORM PATIENT CALCULATIONS  #
##################################

PATIENT_RESULTS = {'PATIENT_ID': [],
                'fSNP_COUNT': [],
                'pNS_fSNP': [],
                'NEUTRAL_pNS_fSNP': [],
                'TRANSITION_STABLE': [],
                'TRANSITION_LOSS': [],
                'TRANSITION_GAIN': [],
                'TRANSITION_ABSENT': [],
                'fSNP_COUNT': [],
                'NON_EFFICACIOUS': []
                }



for ind,_patient in enumerate(PATIENT_DATA['PATIENTS']):
    _data = PATIENT_DATA['UNFIXED_ARRAYS'][ind]
    _absent, _stable, _loss, _gain = \
    serf.transition_matrix_PA(_data, ANNOTATION_DICT, exclude=PATIENT_DATA['MISSING'][ind])

    _effective = int(_patient in ['Patient09', 'Patient10', 'Patient11', 'Patient12'])

    PATIENT_RESULTS['PATIENT_ID'].append(_patient)
    PATIENT_RESULTS['TRANSITION_STABLE'].append(_stable)
    PATIENT_RESULTS['TRANSITION_LOSS'].append(_loss)
    PATIENT_RESULTS['TRANSITION_GAIN'].append(_gain)
    PATIENT_RESULTS['TRANSITION_ABSENT'].append(_absent)
    PATIENT_RESULTS['NON_EFFICACIOUS'].append(_effective)
    PATIENT_RESULTS['fSNP_COUNT'].append(_effective)
    PATIENT_RESULTS['pNS_fSNP'].append(_effective)
    PATIENT_RESULTS['NEUTRAL_pNS_fSNP'].append(_effective)

PATIENT_DF = pd.DataFrame(PATIENT_RESULTS)

PATINET_DF.to_csv('{}{}PATIENT_data.csv'.format(PROCESSED_PATH, os.sep))
