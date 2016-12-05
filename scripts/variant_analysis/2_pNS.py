#!/usr/bin/env python

#Calculate pNS for vSNPs detected at each time point and for each patient.

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
