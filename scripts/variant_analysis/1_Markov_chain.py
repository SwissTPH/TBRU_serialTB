#!/usr/bin/env python

#Use bootstrapping to calculate the transition frequencies for detected vSNPs

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

#############################
#  CALCUALTE MC TRANSITIONS #
#############################

MC_output = {'TRANSITION': [],
             'SNP_TYPE': [],
             'MEDIAN': [],
             '95CI_LOWER': [],
             '95CI_UPPER': [],
             'NON_EFFICACIOUS': [],
             'COUNT': []
            }

transf = lambda x,pos,ind,dn: float(x[pos][ind])/np.sum([x[pos][y] for y in dn])



for patient_type,patient_set in {0:[0,1,2,3,5,6,8,9],1:[4,7,10,11]}.items():
    print('Combined data, inefficacious treatment:{}'.format(patient_type))
    ###COMBINED###
    trans_NS = np.concatenate(tuple(serf.MC_PA_CI(PATIENT_DATA['UNFIXED_ARRAYS'][ind],
                                             ANNOTATION_DICT,
                                             exclude=PATIENT_DATA['MISSING'][ind]) for ind in patient_set))
    MC_PA_trials = np.array([[Counter(sku.resample(trans_NS, random_state=trial)).get(x,0) for x in [0,1,2,-1]]\
                             for trial in range(1000)])
    MC_PA_numerators = np.mean(MC_PA_trials, axis=0)
    MC_PA_denominators = np.mean(MC_PA_trials, axis=0) + np.array([np.mean(MC_PA_trials, axis=0)[y] for y in [3,2,1,0]])
    MC_PA_means = MC_PA_numerators/MC_PA_denominators
    PA_transition_labels = ['0->0','1->1', '1->0', '0->1']

    dns = [[0,3],[1,2],[1,2],[0,3]]
    _tmp_rsl = []
    for ind in [1,2,-1,0]:
        _prs = [transf(MC_PA_trials,x,ind,dns[ind]) for x in range(1000)]
        CI_L, CI_H = np.argsort(_prs)[[24,975]] #95% CI indices
        MC_output['TRANSITION'].append(PA_transition_labels[ind])
        MC_output['SNP_TYPE'].append('COMBINED')
        MC_output['MEDIAN'].append(MC_PA_means[ind])
        MC_output['95CI_LOWER'].append(_prs[CI_L])
        MC_output['95CI_UPPER'].append(_prs[CI_H])
        MC_output['COUNT'].append(Counter(trans_NS).get(ind,0))
        MC_output['NON_EFFICACIOUS'].append(patient_type)
        print('Pr({0}) is {1:4.3f} ({2:4.3f}-{3:4.3f})'.format(PA_transition_labels[ind],MC_PA_means[ind], _prs[CI_L], _prs[CI_H]))


    ###SYNONYMOUS###
    print('Synonymous data, inefficacious treatment:{}'.format(patient_type))
    trans_NS = np.concatenate(tuple(serf.MC_PA_CI(PATIENT_DATA['UNFIXED_ARRAYS'][ind],
                                             ANNOTATION_DICT,
                                             exclude=PATIENT_DATA['MISSING'][ind],
                                             mutation_type='synonymous') for ind in patient_set))
    MC_PA_trials = np.array([[Counter(sku.resample(trans_NS, random_state=trial)).get(x,0) for x in [0,1,2,-1]]\
                             for trial in range(1000)])
    MC_PA_numerators = np.mean(MC_PA_trials, axis=0)
    MC_PA_denominators = np.mean(MC_PA_trials, axis=0) + np.array([np.mean(MC_PA_trials, axis=0)[y] for y in [3,2,1,0]])
    MC_PA_means = MC_PA_numerators/MC_PA_denominators
    PA_transition_labels = ['0->0','1->1', '1->0', '0->1']

    dns = [[0,3],[1,2],[1,2],[0,3]]
    _tmp_rsl = []
    for ind in [1,2,-1,0]:
        _prs = [transf(MC_PA_trials,x,ind,dns[ind]) for x in range(1000)]
        CI_L, CI_H = np.argsort(_prs)[[24,975]] #95% CI indices
        MC_output['TRANSITION'].append(PA_transition_labels[ind])
        MC_output['SNP_TYPE'].append('SYN')
        MC_output['MEDIAN'].append(MC_PA_means[ind])
        MC_output['95CI_LOWER'].append(_prs[CI_L])
        MC_output['95CI_UPPER'].append(_prs[CI_H])
        MC_output['COUNT'].append(Counter(trans_NS).get(ind,0))
        MC_output['NON_EFFICACIOUS'].append(patient_type)
        print('Pr({0}) is {1:4.3f} ({2:4.3f}-{3:4.3f})'.format(PA_transition_labels[ind],MC_PA_means[ind], _prs[CI_L], _prs[CI_H]))

    ###NON-SYNONYMOUS###
    print('Non-synonymous data, inefficacious treatment:{}'.format(patient_type))
    trans_NS = np.concatenate(tuple(serf.MC_PA_CI(PATIENT_DATA['UNFIXED_ARRAYS'][ind],
                                             ANNOTATION_DICT,
                                             exclude=PATIENT_DATA['MISSING'][ind],
                                             mutation_type='nonsynonymous') for ind in patient_set))
    MC_PA_trials = np.array([[Counter(sku.resample(trans_NS, random_state=trial)).get(x,0) for x in [0,1,2,-1]]\
                             for trial in range(1000)])
    MC_PA_numerators = np.mean(MC_PA_trials, axis=0)
    MC_PA_denominators = np.mean(MC_PA_trials, axis=0) + np.array([np.mean(MC_PA_trials, axis=0)[y] for y in [3,2,1,0]])
    MC_PA_means = MC_PA_numerators/MC_PA_denominators
    PA_transition_labels = ['0->0','1->1', '1->0', '0->1']

    dns = [[0,3],[1,2],[1,2],[0,3]]
    _tmp_rsl = []
    for ind in [1,2,-1,0]:
        _prs = [transf(MC_PA_trials,x,ind,dns[ind]) for x in range(1000)]
        CI_L, CI_H = np.argsort(_prs)[[24,975]] #95% CI indices
        MC_output['TRANSITION'].append(PA_transition_labels[ind])
        MC_output['SNP_TYPE'].append('NSY')
        MC_output['MEDIAN'].append(MC_PA_means[ind])
        MC_output['95CI_LOWER'].append(_prs[CI_L])
        MC_output['95CI_UPPER'].append(_prs[CI_H])
        MC_output['COUNT'].append(Counter(trans_NS).get(ind,0))
        MC_output['NON_EFFICACIOUS'].append(patient_type)
        print('Pr({0}) is {1:4.3f} ({2:4.3f}-{3:4.3f})'.format(PA_transition_labels[ind],MC_PA_means[ind], _prs[CI_L], _prs[CI_H]))

MC_df = pd.DataFrame(MC_output)

MC_df.to_csv('{}{}4_MC_patients.csv'.format(PROCESSED_PATH, os.sep))
