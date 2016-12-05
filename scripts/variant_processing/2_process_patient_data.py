#!/usr/bin/env python

#Pre-process sequencing data for serial isolate patients

#NB written for Python 3.5+

################
# DEPENDENCIES #
################

import os
import pickle
import numpy as np
import serial_functions as serf #custom functions (./serial_functions.py)

##########################
#ESTABLISH REFERENCE INFO#
##########################

CWD = os.getcwd() #get current working directory

REFERENCE_PATH = os.path.abspath(os.path.join(CWD, os.pardir, os.pardir, 'data', 'external'))
INTERIM_PATH = os.path.abspath(os.path.join(CWD, os.pardir, os.pardir, 'data', 'interim'))
DATA_PATH = os.path.abspath(os.path.join(CWD, os.pardir, os.pardir, 'data', 'raw'))
FIXED_PATH = os.path.abspath(os.path.join(CWD, os.pardir, os.pardir, 'data', 'raw', 'FIXED_SNP'))
ANNOTATED_PATH = os.path.abspath(os.path.join(CWD, os.pardir, os.pardir, 'data', 'raw', 'ANNOTATED_vSNP'))

REF_ANNOTATIONS = {}
for line in open('{}{}H37Rv_annotation2sytems.ptt'.format(REFERENCE_PATH, os.sep)):
    split = line.strip().split('\t')
    if len(split)<=6 and len(split)>3 and split[0][0].isdigit():
        REF_ANNOTATIONS[int(split[1])-1] = [int(split[2]), split[3], split[4]]
    if len(split)>6 and split[0][0].isdigit():
        REF_ANNOTATIONS[int(split[1])-1] = [int(split[2]), split[3], split[7]]

CCDC_ANNOTATIONS = {}
for line in open('{}{}CCDC5079_annotation2systems.ptt'.format(REFERENCE_PATH, os.sep)):
    split = line.strip().split('\t')
    if len(split)<=6 and len(split)>3 and split[0][0].isdigit():
        CCDC_ANNOTATIONS[int(split[1])-1] = [int(split[2]), split[3], split[4]]
    if len(split)>6 and split[0][0].isdigit():
        CCDC_ANNOTATIONS[int(split[1])-1] = [int(split[2]), split[3], split[7]]


LOCUS_MATRIX = np.array([[k,v[0]] for (k,v) in REF_ANNOTATIONS.items()], dtype=float)
LOCUS_MATRIX = LOCUS_MATRIX[np.argsort(LOCUS_MATRIX[:,0])]

LOCUS_MATRIX_C = np.array([[k,v[0]] for (k,v) in CCDC_ANNOTATIONS.items()], dtype=float)
LOCUS_MATRIX_C = LOCUS_MATRIX_C[np.argsort(LOCUS_MATRIX_C[:,0])]

LOCI = np.array([REF_ANNOTATIONS[int(x)][-1] for x in LOCUS_MATRIX[:,0]])
SIZE_DICT = {x:float(LOCUS_MATRIX[ind][1])-LOCUS_MATRIX[ind][0] for (ind,x) in enumerate(LOCI)}
MTB_GENOME = ''.join([line.strip() for line in open('{}{}MTB_anc.fa'.format(REFERENCE_PATH, os.sep)) if line[0]!='>'])
CCDC5079_GENOME = ''.join([line.strip() for line in open('{}{}NC_021251.fa'.format(REFERENCE_PATH, os.sep)) if line[0]!='>'])

CCDC5079_SIZE = len(CCDC5079_GENOME) #based on NC_021251 (CCDC5079)
CCDC_codons, CCDC_total_codons, CCDC_total_features = serf.codon_counter('{}{}NC_021251_annotations.fasta'.format(REFERENCE_PATH, os.sep))
CCDC_codon_pool = serf.generate_codon_pool(CCDC_codons, CCDC_total_features)


CODON_TABLE = {'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAT': 'N', 'ACA': 'T',
        'ACC': 'T', 'ACG': 'T', 'ACT': 'T', 'AGA': 'R', 'AGC': 'S', 'AGG': 'R',
        'AGT': 'S', 'ATA': 'I', 'ATC': 'I', 'ATG': 'M', 'ATT': 'I', 'CAA': 'Q',
        'CAC': 'H', 'CAG': 'Q', 'CAT': 'H', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P',
        'CCT': 'P', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R', 'CTA': 'L',
        'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 'GAA': 'E', 'GAC': 'D', 'GAG': 'E',
        'GAT': 'D', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A', 'GGA': 'G',
        'GGC': 'G', 'GGG': 'G', 'GGT': 'G', 'GTA': 'V', 'GTC': 'V', 'GTG': 'V',
        'GTT': 'V', 'TAA': 'STOP', 'TAC': 'Y', 'TAG': 'STOP', 'TAT': 'Y',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'TGA': 'STOP', 'TGC': 'C',
        'TGG': 'W', 'TGT': 'C', 'TTA': 'L', 'TTC': 'F', 'TTG': 'L', 'TTT': 'F'}

#######################
#PARSE SEQUENCING DATA#
#######################

DATA_DICT = {} #recipient dictionary
ANNOTATION_DICT = {} # ANNOTATION_DICT['4368405'] > [['Rv3885c', 'synonymous', '39', 'Val/V-Val/V', 'GTC-GTT']]
FIXED_ANNOTATION_DICT = {} # FIXED_ANNOTATION_DICT['1532463'] > [['CFBS07230', 'synonymous', 240, 'G/G', 'GGG-GGT', 0]]
FIXED_MUTATION_TYPE={}

for line in open('{}{}SAMPLES.list'.format(DATA_PATH, os.sep)):
    if line[0] == '>':
        _patient = line.strip()[1:]
        DATA_DICT[_patient]=[]
        FIXED_MUTATION_TYPE[_patient]=[]
    if line[2] == '.':
        line=line.strip()
        _posfreqs = []
        _posclass = []
        for entry in open('{}{}'.format(FIXED_PATH,line)):
            _data = entry.split()
            if _data[0][0]=='%': _data[0]=_data[0][1:]
            _posfreqs.append(_data[0])
            _posclass.append(serf.mutation_classifier(_data[1], _data[2]))
            if _data[0] in FIXED_ANNOTATION_DICT:
            #    FIXED_ANNOTATION_DICT[_data[0]].append(serf.SNP_annotate(int(_data[0])-1, _data[1], _data[2], locus_matrix=LOCUS_MATRIX, reference=REF_ANNOTATIONS, genome=MTB_GENOME, codon_table=CODON_TABLE))
                FIXED_ANNOTATION_DICT[_data[0]].append(serf.SNP_annotate(int(_data[0])-1, _data[1], _data[2], locus_matrix=LOCUS_MATRIX_C, reference=CCDC_ANNOTATIONS, genome=CCDC5079_GENOME, codon_table=CODON_TABLE))
            if _data[0] not in FIXED_ANNOTATION_DICT:
                FIXED_ANNOTATION_DICT[_data[0]]=[serf.SNP_annotate(int(_data[0])-1, _data[1], _data[2], locus_matrix=LOCUS_MATRIX_C, reference=CCDC_ANNOTATIONS, genome=CCDC5079_GENOME, codon_table=CODON_TABLE)]
        DATA_DICT[_patient].append(['fixed']+_posfreqs)
        FIXED_MUTATION_TYPE[_patient]+=_posclass
    if line[2].isdigit():
        line=line.strip()
        _posfreqs = []
        for entry in open('{}{}{}'.format(ANNOTATED_PATH,os.sep,line)):
            _data = entry.split()
            #if float(_data[3][:-1])>1.5: #Filter SNPs with a frequency of less than 1.5%
            if _data[0][0]=='%': _data[0]=_data[0][1:]
            _posfreqs.append((_data[0],float(_data[3][:-1])/100))
            if _data[0] in ANNOTATION_DICT:
                MUTATION_CLASS = serf.mutation_classifier(_data[1], _data[2]) #classify mutations as transitions (0) of transversions (1)
                ANNOTATION_DICT[_data[0]].append([_data[8]]+[_data[5]]+[_data[4]]+[_data[6]]+[_data[7]]+[MUTATION_CLASS])
            if _data[0] not in ANNOTATION_DICT:
                MUTATION_CLASS = serf.mutation_classifier(_data[1], _data[2]) #classify mutations as transitions (0) of transversions (1)
                ANNOTATION_DICT[_data[0]]=[[_data[8]]+[_data[5]]+[_data[4]]+[_data[6]]+[_data[7]]+[MUTATION_CLASS]]
        DATA_DICT[_patient].append([line[2:4]]+_posfreqs)


FIXED_ARRAY = np.zeros((2,13)) #Need to start with a 2x13 array to be able to iterate correctly.

for ind,k in enumerate(sorted(DATA_DICT.keys())):
    for _locus in DATA_DICT[k][0][1:]:
        if float(_locus) in FIXED_ARRAY[:,0]:
            _mapper = list(FIXED_ARRAY[:,0])
            _bob = _mapper.index(float(_locus))
            FIXED_ARRAY[_bob,ind+1]+=1
        if _locus not in FIXED_ARRAY[:,0]:
            _new = np.zeros(13)
            _new[0]+=int(_locus)
            _new[ind+1]+=1
            FIXED_ARRAY = np.vstack((FIXED_ARRAY, _new))

#delete the first two dummy lines:
FIXED_ARRAY=FIXED_ARRAY[2:]

#Get the data for the re-samplings of Patient12 t=0:

_values = DATA_DICT['Patient12'][1:5]

p12_0 = np.zeros((2,5))

for _ind,_v in enumerate(_values):
        for (_locus,_freq) in _v[1:]:
            if float(_locus) in p12_0[:,0]:
                _mapper = list(p12_0[:,0])
                _bob = _mapper.index(float(_locus))
                p12_0[_bob,_ind+1]+=float(_freq)
            if float(_locus) not in p12_0[:,0]:
                _new = np.zeros(5)
                _new[0]+=int(_locus)
                _new[_ind+1]+=float(_freq)
                p12_0 = np.vstack((p12_0, _new))

p12_0 = p12_0[2:] #delete top two dummy lines

##########################
#GLOBAL PATIENT VARIABLES#
##########################

p1 = serf.populate_poly(DATA_DICT, patient='Patient01')
p2 = serf.populate_poly(DATA_DICT, patient='Patient02')
p3 = serf.populate_poly(DATA_DICT, patient='Patient03')
p4 = serf.populate_poly(DATA_DICT, patient='Patient04')
p5 = serf.populate_poly(DATA_DICT, patient='Patient05')
p6 = serf.populate_poly(DATA_DICT, patient='Patient06')
p7 = serf.populate_poly(DATA_DICT, patient='Patient07')
p8 = serf.populate_poly(DATA_DICT, patient='Patient08')
p9 = serf.populate_poly(DATA_DICT, patient='Patient09')
p10 = serf.populate_poly(DATA_DICT, patient='Patient10')
p11 = serf.populate_poly(DATA_DICT, patient='Patient11')
p12 = serf.populate_poly(DATA_DICT, patient='Patient12')

PATIENT_DATA = {}
PATIENT_DATA['UNFIXED_ARRAYS'] = [
                p1, p2, p6, p7, p11, p3,
                p4, p12, p5, p8, p9, p10
                ]

PATIENT_DATA['MGIT'] = [6, 2, 8, 8, 4, 8, 4, 20, 6, 26, 16, 16]

PATIENT_DATA['DS'] = ['g','g','b','b','r','g','g','k','g','b','r','r']

PATIENT_DATA['DS2'] = ['g','g','b','b','k','g','g','k','g','b','r','r']

PATIENT_DATA['TP'] = [
                'DS','DS','Hr','Hr','MDR-re',
                'DS','DS','MDR-re','DS','Hr',
                'MDR','MDR'
                ]

PATIENT_DATA['TP2'] = [
                'DS','DS','Hr','Hr','MDR',
                'DS','DS','MDR','DS','Hr',
                'MDR','MDR'
                ]

PATIENT_DATA['TP_SORTING'] = np.array([0,1,5,6,8,2,3,9,10,11,4,7])

PATIENT_DATA['PATIENTS'] = [
                'Patient01', 'Patient02', 'Patient06','Patient07',
                'Patient11', 'Patient03','Patient04','Patient12',
                'Patient05', 'Patient08','Patient09', 'Patient10'
                ] #sorted(DATA_DICT.keys())

PATIENT_DATA['TIMEPOINTS'] = ['00','02','04','06','08','16','24']

PATIENT_DATA['TM'] = [0,2,4,6,8,16,24]

PATIENT_DATA['MISSING'] = [
                None, None, None, None,
                None, 3, 2, None, None,
                None, 5, 2
                ]
PATIENT_DATA['DEPTH'] = {'Patient01': {0:1106.8, 2:1482.9, 4:1266.7, 6:1515.6, 8:2552.0},
                         'Patient02': {0:1424.5, 2:1360.2, 4:1269.4, 6:1454.0},
                         'Patient03': {0:1311.7, 2:1039.1, 6: 708.8, 8:1284.8},
                         'Patient04': {0: 948.0, 4:1097.8},
                         'Patient05': {0:1026.5, 2:1059.2, 4:1018.3, 6: 956.4, 8:1562.6},
                         'Patient06': {0:1388.3, 2: 353.4, 4:1431.0, 6:1339.4, 8:1308.3},
                         'Patient07': {0:1398.2, 2:1093.7, 4:1683.1, 6: 855.1, 8: 991.8},
                         'Patient08': {0: 775.9, 2:1350.7, 4:1425.0, 6:1356.9, 8:1002.2},
                         'Patient09': {0:1237.1, 2: 858.0, 4:1532.2, 6:1718.9},
                         'Patient10': {0:2101.0, 4:1262.8, 6: 940.8, 8:1372.6},
                         'Patient11': {0:1664.6, 2: 700.2, 4:1485.8, 6:1431.6, 8:1268.8},
                         'Patient12': {0:1724.7, 2: 904.3, 4: 822.6, 6:1790.7, 8:1422.4}
                         }

###############################
#DUMP KEY VARIABLES TO PICKLES#
###############################

pickle.dump(PATIENT_DATA, open('{}{}2_patient_data.pkl'.format(INTERIM_PATH, os.sep),'wb'))
pickle.dump(ANNOTATION_DICT, open('{}{}2_unfixed_annotation.pkl'.format(INTERIM_PATH, os.sep),'wb'))
pickle.dump(FIXED_ANNOTATION_DICT, open('{}{}2_fixed_annotation.pkl'.format(INTERIM_PATH, os.sep),'wb'))
pickle.dump(DATA_DICT, open('{}{}2_unfixed_mutations.pkl'.format(INTERIM_PATH, os.sep),'wb'))
pickle.dump(FIXED_ARRAY, open('{}{}2_fixed_mutations.pkl'.format(INTERIM_PATH, os.sep),'wb'))
