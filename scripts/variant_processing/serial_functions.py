#!/usr/bin/env python

#Custom analysis functions for processing data from Serial Isolates of TB patients

#NB written for Python 3.5+

################
# DEPENDENCIES #
################

import numpy as np
import sklearn.utils as sku
import scipy.stats as ss
from collections import Counter
import random

#############
# FUNCTIONS #
#############

def populate_poly(data_dict, patient='Patient02'):
    """
    Extract SNP frequencies from raw data.

    INPUTS:
    -------
    patient: str, dict key for data_dict
    data_dict: dict, {PatientID: [['timepoint', (locus, frequency)]...]}

    OUTPUT:
    -------
    nd.array, first column are loci, sequenital columns are frequencies
    """
    _polys = np.zeros((2,8)) #Need to start with a 2x8 array to be able to iterate correctly.

    _timepoints = ['00','02','04','06','08','16','24']

    _values = data_dict[patient][1:]

    for _v in _values:
        for (_locus,_freq) in _v[1:]:
            _ind = _timepoints.index(_v[0])
            if float(_locus) in _polys[:,0]:
                _mapper = list(_polys[:,0])
                _bob = _mapper.index(float(_locus))
                _polys[_bob,_ind+1]+=float(_freq)
            if float(_locus) not in _polys[:,0]:
                _new = np.zeros(8)
                _new[0]+=int(_locus)
                _new[_ind+1]+=float(_freq)
                _polys = np.vstack((_polys, _new))

    _polys = _polys[2:] #Remove first two rows.
    return _polys

def polymorphic_sites_plot(data, cutoff=0.005, colour='red', return_data=False, savefig=False):
    """Plot of v-SNP counts.

    INPUT:
    ------
    data: 2D-array generated from populate_poly()
    cutoff: float, minimal frequency required for counting
    colour: color of line, can take any value accepted by matplotlib.
    return_data: Boolean, return the plotted table
    savefig: give string with filename to generate if you want to save the figure.
    """
    _timepoints = np.array([0,2,4,6,8,16,24])

    _polymorphic_sites = np.array([len(np.where(data[:,x]>cutoff)[0]) for x in np.arange(1,8)])
    #only plot (consider) timepoints for which we have data.
    _to_consider = np.nonzero(np.array([len(np.where(data[:,x]>0)[0]) for x in np.arange(1,8)]))
    _x, _y = _timepoints[_to_consider], _polymorphic_sites[_to_consider]
    c_label = 'Cutoff = %.1f' %(cutoff*100)+'%'
    plt.plot(_x, _y, color=colour, label=c_label,lw=1.5)

    if savefig!=False:
        plt.savefig(savefig)

    if return_data!=False:
        return _x, _y

def SFS(data, title='Patient__', plots=[0,2,4,6], grid=(2,2)):
     _timepoints = [0,2,4,6,8,16,24]

     for _sn,_pl in enumerate(plots):
        _i = _timepoints.index(_pl)
        if grid!=False:
            plt.subplot(grid[0],grid[1],_sn+1)
            plt.title('Week %s' %_timepoints[_i])
        _to_plot = data[:,_i+1][data[:,_i+1]>0]
        plt.hist(_to_plot,bins=np.arange(0,1.01,0.01), label='Week %s' %_timepoints[_i], alpha=0.5)
        plt.xlim(0,1)

def SNP_fate(data, cutoff=0.005, method='imshow'):
    _a = np.zeros(np.shape(data[:,1:]))
    _a[data[:,1:]>cutoff]+=1
    _relevant = np.where(_a.sum(axis=1)>1)
    _labs = data[:,0][_relevant]
    _tm = np.array([0,2,4,6,8,16,24])
    if method=='imshow':
        plt.imshow(_a[_relevant], interpolation='none',cmap='Greys')
        plt.yticks(np.arange(1,len(_relevant)),_labs)
    if method=='plot':
        if len(_relevant)>0:
            for _ind,_ys in enumerate(data[:,1:][_relevant]):
                plt.plot(_tm[_a.sum(axis=0)>0], _ys[_a.sum(axis=0)>0], label='SNP at %i' %data[:,0][_relevant][_ind],lw=1.5)

def SNP_fate_NSI(data, annotation_dict, cutoff=0.005, method='imshow'):
    _a = np.zeros(np.shape(data[:,1:]))
    _a[data[:,1:]>cutoff]+=1
    _relevant = np.where(_a.sum(axis=1)>1)
    _labs = data[:,0][_relevant]
    _tm = np.array([0,2,4,6,8,16,24])
    if method=='imshow':
        plt.imshow(_a[_relevant], interpolation='none',cmap='Greys')
        plt.yticks(np.arange(1,len(_relevant)),_labs)
    if method=='plot':
        if len(_relevant)>0:
            for _ind,_ys in enumerate(data[:,1:][_relevant]):
                _mutation_type_color = {'synonymous':'grey', 'nonsynonymous':'red','IGR':'dodgerblue'}
                _mutation_type = _mutation_type_color.get(annotation_dict[str(int(data[:,0][_relevant][_ind]))][0][1],'yellow')
                plt.plot(_tm[_a.sum(axis=0)>0], _ys[_a.sum(axis=0)>0],color=_mutation_type, label='SNP at %i' %data[:,0][_relevant][_ind],lw=1.5,alpha=0.8)

def SNP_fate_DR(data, annotation_dict, cutoff=0.005, method='imshow'):
    _a = np.zeros(np.shape(data[:,1:]))
    _a[data[:,1:]>cutoff]+=1
    _relevant = np.where(_a.sum(axis=1)>1)
    _labs = data[:,0][_relevant]
    _tm = np.array([0,2,4,6,8,16,24])
    if method=='imshow':
        plt.imshow(_a[_relevant], interpolation='none',cmap='Greys')
        plt.yticks(np.arange(1,len(_relevant)),_labs)
    if method=='plot':
        if len(_relevant)>0:
            for _ind,_ys in enumerate(data[:,1:][_relevant]):
                _mutation_type_color = ['grey','red']
                _mutation_type = _mutation_type_color[int(annotation_dict[str(int(data[:,0][_relevant][_ind]))][0][0] in DR_set)]
                _mutation_type_style = {'nonsynonymous':'solid', 'synonymous':'dashed','IGR':'dotted'}
                _mutation_style = _mutation_type_style.get(annotation_dict[str(int(data[:,0][_relevant][_ind]))][0][1],'dashdot')
                plt.plot(_tm[_a.sum(axis=0)>0], _ys[_a.sum(axis=0)>0],color=_mutation_type, linestyle=_mutation_style,label='SNP at %i' %data[:,0][_relevant][_ind],lw=1.5,alpha=0.6)


def confidence_interval(data, CI=0.95, normal=False):
        """Calculate the confidence interval for a dataset
        which is distribution blind.

        Can be used in the case of a normal distribution expectation.

        INPUT:
        ------
        data: 1D np.array
        CI: float < 1.0, confidence interval of interest
            expects 0.95 for 95% CI, 0.9 for 90% CI...
        normal: bool, determine CI based on normal distribution.


        OUTPUT:
        -------
        lower_CI, upper_CI: float


        NOTES:
        ------
        Designed to be used in conjuncion with bootstrapping

        Adapted from:
        http://stackoverflow.com/questions/19124239/scikit-learn-roc-curve-with-confidence-intervals
        http://adventuresinpython.blogspot.ch/2012/12/confidence-intervals-in-python.html

        """

        _CI = (1.+CI)/2.

        _lower_CI = sorted(data)[int((1.-_CI)*len(data))]
        _upper_CI = sorted(data)[int((_CI)*len(data))]

        _R = (_lower_CI, _upper_CI)

        if normal!=False:
            _n, _min_max, _mean, _var, _skew, _kurt = ss.describe(data)
            _std=np.sqrt(_var)

            _R = ss.norm.interval(1-_CI,loc=_mean,scale=_std)

        return _R

def simpson_index(data):
    """Calculate Simpson's index from data.

    In ecology this index corresponds to the
    probabilty that two individuals drawn at
    random from a population belong to the
    same group.

    Index takes a value between 0 and 1. The
    greater the value, the more homogeneous a
    sample is.

    1-SI = Simpson's index of diversity (SID)

        Probability of two individuals belonging
        to different groups.

    1/SI = Simpson's reciprocal index

        Effective number of groups in a population.
        This parameters takes a value between 1 and
        K where K is the actual number of groups in
        the sample.

    INPUT:
    ------
    data: 1d array (sums to 1.)

    OUTPUT:
    -------
    Simpason's index: float

    NOTES:
    ------
    Based on http://en.wikipedia.org/wiki/Diversity_index and
    http://www.countrysideinfo.co.uk/simpsons.htm.
    """

    return np.sum(data**2)

def euclidean_distance(data1, data2):
    """Calculate the euclidian distance between two points.

    INPUTS:
    -------
    data1, data2: 1D-array of equal length

    OUTPUT:
    -------
    float

    NOTES:
    ------
    Based on http://en.wikipedia.org/wiki/Euclidean_distance
    """

    return np.sqrt(np.sum((np.array(data1)-np.array(data2))**2))

def heterozygosity(het_calls, genome_size):
    """Calcualtes the heterozygosity of a sample from
    the relative frequencies of the heterozygous alleles
    in the population.

    In this iteration it assumes two alleles per position.

    So Fa = frequency of allele 1, while 1-Fa is frequency
    of allele 2.

    INPUTS:
    -------
    het_calls: 1-D array, values range 0-1.
    genome_size: float, size of genome in question.

    OUTPUT:
    -------
    float

    NOTES:
    ------
    Based on the method for determining heterozygosity in
    Cuevas et al, MBE 2015
    http://www.ncbi.nlm.nih.gov/pubmed/25660377
    """

    _N = len(het_calls)

    _het1 = het_calls
    _het2 = np.ones(_N)-het_calls

    _hets = np.sum(np.ones(_N)-(_het1**2+_het2**2))

    return _hets/float(genome_size)


def codon_counter(annotation_file):
    """Calcuates codon frequencies in
    an annotated genome.

    Assumes a file where:

    >gene_annotation
    TTGACCGATGACCCCGGTTC...

    the coding sequence is given in a
    single line.

    ARGUMENT:
    ---------
    Text file formatted as described above

    OUTPUT:
    -------
    codon_freq: codon counts, dict
    total_codons: absolute codon count, float
    total_features: absolute feature count, float

    """

    _codon_dict = {}
    _codon_clicker = 0.
    _feature_clicker = 0.

    for _line in open(annotation_file):
        if _line[0]!='>':
            _feature_clicker+=1.
            _sequence = _line.strip()
            for i in range(0,len(_sequence),3):
                _codon = _sequence[i:i+3]
                _codon_clicker+=1.
                if _codon in _codon_dict:
                    _codon_dict[_codon]+=1
                if _codon not in _codon_dict:
                    _codon_dict[_codon]=1

    return _codon_dict, _codon_clicker, _feature_clicker

def generate_codon_pool(codon_dict, N):
    """Generate a weighted pool of codons to
    be used for codon resampling from a fake
    genome-defined average gene.

    ARGUMENTS:
    ----------
    codon_dict: dictionary of codon counts in
                a genome from codon_counter()

    N: number of genes in the genome.
        total_features from codon_counter()

    OUTPUT:
    -------
    codon_pool: list

    NOTES:
    ------
    Assumes numpy is imported

    """

    _codon_pool = []

    for _key in sorted(codon_dict.keys()):
        weight = int(np.ceil(codon_dict[_key]/N))
        _codon_pool+=[_key,]*weight

    return _codon_pool

def codon_mutate(codon):
    """Mutates codon at one position.
    Checks whether the mutations was
    synonymous or nonsynonymous.

    ARGUMENTS:
    ----------
    codon: str (3 bases)

    OUTPUT:
    -------
    0,1 for S or NS respectively
    """

    codon_table = {'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAU': 'N', 'ACA': 'T',
        'ACC': 'T', 'ACG': 'T', 'ACU': 'T', 'AGA': 'R', 'AGC': 'S', 'AGG': 'R',
        'AGU': 'S', 'AUA': 'I', 'AUC': 'I', 'AUG': 'M', 'AUU': 'I', 'CAA': 'Q',
        'CAC': 'H', 'CAG': 'Q', 'CAU': 'H', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P',
        'CCU': 'P', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R', 'CUA': 'L',
        'CUC': 'L', 'CUG': 'L', 'CUU': 'L', 'GAA': 'E', 'GAC': 'D', 'GAG': 'E',
        'GAU': 'D', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A', 'GGA': 'G',
        'GGC': 'G', 'GGG': 'G', 'GGU': 'G', 'GUA': 'V', 'GUC': 'V', 'GUG': 'V',
        'GUU': 'V', 'UAA': 'STOP', 'UAC': 'Y', 'UAG': 'STOP', 'UAU': 'Y',
        'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S', 'UGA': 'STOP', 'UGC': 'C',
        'UGG': 'W', 'UGU': 'C', 'UUA': 'L', 'UUC': 'F', 'UUG': 'L', 'UUU': 'F'}

    #Format the sequnce.
    codon = codon.upper()

    if 'T' in codon:
        codon = codon.replace('T', 'U')

    _bases = ['A', 'C', 'G', 'U'] #all bases

    _mut_pos = np.random.randint(3) #pick codon position to mutate

    _bases.remove(codon[_mut_pos]) #remove the ancestral base

    _mut_base = _bases[np.random.randint(3)] #pick random base

    _alt_codon = codon[:_mut_pos]+_mut_base+codon[_mut_pos+1:] #mutate codon

    if codon_table[codon]==codon_table[_alt_codon]: return 0
    if codon_table[codon]!=codon_table[_alt_codon]: return 1

def codon_mutate_TiTv_GC(codon, gc_content=0.5, ti=0.5):
    """Mutates codon at one position.
    Checks whether the mutations was
    synonymous or nonsynonymous.

    ARGUMENTS:
    ----------
    codon: str (3 bases)
    gc_content: float, if specified will be used for
                for calculating the substitution_matrix.
                otherwise will assume 0.5.
    ti: float, transition probability

    OUTPUT:
    -------
    0,1 for S or NS respectively

    NOTES:
    ------
    Substitution matrix calculated based on
    Jukes-Cantor 1969 in the basic (default) case.
    Changing 'ti' modifies the approach according to
    Kimura 1980. Changing 'gc_content' modifies it based
    on Felsenstein 1981. Changing both 'ti' and 'gc_content'
    approximates Hasegawa, Kishino and Yano 1985.
    """

    codon_table = {'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAU': 'N', 'ACA': 'T',
        'ACC': 'T', 'ACG': 'T', 'ACU': 'T', 'AGA': 'R', 'AGC': 'S', 'AGG': 'R',
        'AGU': 'S', 'AUA': 'I', 'AUC': 'I', 'AUG': 'M', 'AUU': 'I', 'CAA': 'Q',
        'CAC': 'H', 'CAG': 'Q', 'CAU': 'H', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P',
        'CCU': 'P', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R', 'CUA': 'L',
        'CUC': 'L', 'CUG': 'L', 'CUU': 'L', 'GAA': 'E', 'GAC': 'D', 'GAG': 'E',
        'GAU': 'D', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A', 'GGA': 'G',
        'GGC': 'G', 'GGG': 'G', 'GGU': 'G', 'GUA': 'V', 'GUC': 'V', 'GUG': 'V',
        'GUU': 'V', 'UAA': 'STOP', 'UAC': 'Y', 'UAG': 'STOP', 'UAU': 'Y',
        'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S', 'UGA': 'STOP', 'UGC': 'C',
        'UGG': 'W', 'UGU': 'C', 'UUA': 'L', 'UUC': 'F', 'UUG': 'L', 'UUU': 'F'}

    #Format the sequnce.
    codon = codon.upper()

    if 'T' in codon:
        codon = codon.replace('T', 'U')

    _pa = _pt = (1-gc_content)/2. #genomewide probability of base being A/T
    _pg = _pc = gc_content/2.  #genomewide probability of base being C/G

    _a = ti #null expectation for transition
    _b = 1-ti #null expectation for transversion

    _sub_matrix = {
    'A':[_b*_pc/((2*_b+_a)*_pc), _a*_pg/((2*_b+_a)*_pg), _b*_pt/((2*_b+_a)*_pt)],
    'C':[_b*_pa/((2*_b+_a)*_pa), _b*_pg/((2*_b+_a)*_pg), _a*_pt/((2*_b+_a)*_pt)],
    'G':[_a*_pa/((2*_b+_a)*_pa), _b*_pc/((2*_b+_a)*_pc), _b*_pt/((2*_b+_a)*_pt)],
    'U':[_b*_pa/((2*_b+_a)*_pa), _a*_pc/((2*_b+_a)*_pc), _b*_pg/((2*_b+_a)*_pg)]
    }

    _bases = ['A', 'C', 'G', 'U'] #all bases

    _mut_pos = np.random.randint(3) #pick codon position to mutate

    _sub_params = _sub_matrix[codon[_mut_pos]] #Get substitution probabilities
    _bases.remove(codon[_mut_pos]) #remove the ancestral base

    _adjusted_bases = []
    for _ind,_val in enumerate(_bases):
        _expansion = [_val,]*int(_sub_params[_ind]*100)
        _adjusted_bases += _expansion

    _mut_base = _adjusted_bases[np.random.randint(len(_adjusted_bases))] #pick random base

    _alt_codon = codon[:_mut_pos]+_mut_base+codon[_mut_pos+1:] #mutate codon

    if codon_table[codon]==codon_table[_alt_codon]: return 0
    if codon_table[codon]!=codon_table[_alt_codon]: return 1

def calculate_NS(sequence):
    """Calculates the synonymous/nonsynonymous score
    for a sequence.

    It takes the sequence, generates all codons and
    calculates all that are 1 mutation away. Then it
    scores the resulting population for mutational
    potential.

    Used for simplified dN/dS calculations.

    ARGUMENT:
    ---------
    sequence : any nucleotide sequence, at the core should
        be a codon.

    OUTPUT:
    -------
    [(N,S)]
    N : nonsynonymous score per codon
    S : synonymous score per codon

    Note, the function will return a list of tuples if multiple codons are present
    """

    codon_table = {'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAT': 'N', 'ACA': 'T',
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

    #Format the sequnce.
    sequence = sequence.upper()

    if 'T' in sequence:
        sequence = sequence.replace('T', 'U')

    bases = ['A', 'C', 'G', 'U'] #all bases

    sequence_len = len(sequence)
    recipient_list = []

    for codon in [sequence[c:c+3] for c in np.arange(0,len(sequence),3)]: #codonise
        original_aa = codon_table[codon] #expected codon
        codon_N = 0.
        codon_S = 0.
        for ind,base in enumerate(codon):
            codon_marked = codon[:ind]+'X'+codon[ind+1:] #Put 'X' in place of the base to be mutated, this is to avoid replacing all the bases of a certain type in a string...
            alt_codons = [codon_marked.replace('X',b) for b in [x for x in bases if x not in base]] # nested list comprehension, the first gets all other bases, the second generates all codons one mutation apart at given position

            syn = [codon_table[ac] for ac in alt_codons].count(original_aa)/3. #translate alt_codons and see how many match the original aa
            nsyn = 1.-syn #how many dont match

            codon_N+=nsyn
            codon_S+=syn
        recipient_list.append((codon_N, codon_S))

    if len(recipient_list)==1:
        return recipient_list[0][0], recipient_list[0][1]

    if len(recipient_list)>1:
        return recipient_list

def drift_simulation(codon_pool, N_codons, N=1000, method='binomial', pars=False):
    """Simulate a resamplng of a codon pool N times
    each time fishing out N_codons. A binomial random
    draw is used to "mutate" the codon to model ranom
    mutations. The wildtype codons are concatenated
    and the expected number of synonymous and non-
    synonymous mutations is calculated using
    calculate_NS().

    ARGUMENTS:
    ----------
    codon_pool: weighted list of codons done with
                generate_codon_pool()

    N_codons: number of codons to fish out

    N: number of iterations

    method ['derived'|'binomial']: codon mutation model

    pars: return dN and dS separately - in case there are too
          mutations to effectively calculate dN/dS

    OUTPUT:
    -------
    dN/dS values, np.array

    NOTES:
    ------
    This function has many dependencies, ensure that
    all the relevant functions are defined first.

    The probabilites of N/S are defined by _pr_N
    """

    #Define NS probabilities for binomial draws.
    #This is based on the codon table.

    _pr_N = {'AAA': 0.888888888888889, 'AAC': 0.888888888888889, 'AAG': 0.888888888888889, 'AAT': 0.888888888888889,
            'ACA': 0.6666666666666666, 'ACC': 0.6666666666666666, 'ACG': 0.6666666666666666, 'ACT': 0.6666666666666666,
            'AGA': 0.7777777777777778, 'AGC': 0.888888888888889, 'AGG': 0.7777777777777778, 'AGT': 0.888888888888889,
            'ATA': 0.7777777777777778, 'ATC': 0.7777777777777778, 'ATG': 1.0, 'ATT': 0.7777777777777778,
            'CAA': 0.888888888888889, 'CAC': 0.888888888888889, 'CAG': 0.888888888888889, 'CAT': 0.888888888888889,
            'CCA': 0.6666666666666666, 'CCC': 0.6666666666666666, 'CCG': 0.6666666666666666, 'CCT': 0.6666666666666666,
            'CGA': 0.5555555555555556, 'CGC': 0.6666666666666666, 'CGG': 0.5555555555555556, 'CGT': 0.6666666666666666,
            'CTA': 0.5555555555555556, 'CTC': 0.6666666666666666, 'CTG': 0.5555555555555556, 'CTT': 0.6666666666666666,
            'GAA': 0.888888888888889, 'GAC': 0.888888888888889, 'GAG': 0.888888888888889, 'GAT': 0.888888888888889,
            'GCA': 0.6666666666666666, 'GCC': 0.6666666666666666, 'GCG': 0.6666666666666666, 'GCT': 0.6666666666666666,
            'GGA': 0.6666666666666666,'GGC': 0.6666666666666666, 'GGG': 0.6666666666666666, 'GGT': 0.6666666666666666,
            'GTA': 0.6666666666666666, 'GTC': 0.6666666666666666, 'GTG': 0.6666666666666666, 'GTT': 0.6666666666666666,
            'TAA': 0.7777777777777778, 'TAC': 0.888888888888889, 'TAG': 0.888888888888889, 'TAT': 0.888888888888889,
            'TCA': 0.6666666666666666, 'TCC': 0.6666666666666666, 'TCG': 0.6666666666666666, 'TCT': 0.6666666666666666,
            'TGA': 0.888888888888889, 'TGC': 0.888888888888889, 'TGG': 1.0, 'TGT': 0.888888888888889,
            'TTA': 0.7777777777777778, 'TTC': 0.888888888888889, 'TTG': 0.7777777777777778, 'TTT': 0.888888888888889}

    #define recipient variables

    _dNdS = []
    _dNdS_par = []

    for _trial in range(N):#loop over number of trials

        _Ns = 0. #recipient for Nonsynonymous mutations
        _CONCATENATED_SEQUENCE = '' #recipient for concatenated sequence

        for _fish in np.random.randint(0,len(codon_pool),N_codons): #select random codon weighted by genomewide codon usage
            _fished_codon = codon_pool[_fish] #define the fished codon
            _CONCATENATED_SEQUENCE+=_fished_codon #add fished codon to the concatenated sequence
            if method=='derived': _Ns+=codon_mutate(_fished_codon) #use the codon_mutate() to model mutation
            if method=='binomial': _Ns+=np.random.binomial(1, _pr_N[_fished_codon]) #use a binomial random trial to determine if the mutation is Ns

        _eNS = calculate_NS(_CONCATENATED_SEQUENCE)
        _eN = 0.
        _eS = 0.
        for (_a,_b) in _eNS:
            _eN+=_a
            _eS+=_b
        _dNdS.append((_Ns/_eN)/((N_codons-_Ns)/_eS))
        _dNdS_par.append((_Ns/_eN,(N_codons-_Ns)/_eS))

    if pars: return _dNdS_par
    else: return _dNdS

def drift_null_Mtb(sequence, N=1000, method='binomial', pars=False, gc=0.656, ti=0.729):
    """Generate an estimate of the null expectation for
    the pNS of a given Mtb sequence.

    ARGUMENTS:
    ----------
    sequence: str, concatenated codons of interest

    N: number of iterations

    method ['derived'|'binomial']: codon mutation model

    pars: return dN and dS separately - in case there are too
          mutations to effectively calculate dN/dS

    gc: float, gc content of the genome

    ti: float, estimated probability of transitions.

    OUTPUT:
    -------
    expected pNS value under drift, np.array

    NOTES:
    ------
    This function has many dependencies, ensure that
    all the relevant functions are defined first.

    The probabilites of N/S are defined by _pr_N
    """

    #Define NS probabilities for binomial draws.
    #This is based on the codon table.

    _pr_N = {'ACC': 0.66410000000000013, 'ATG': 1.0, 'AAG': 0.80871999999999999,
             'AAA': 0.80523999999999984, 'ATC': 0.73629999999999984, 'AAC': 0.80698000000000003,
             'ATA': 0.85933999999999988, 'AGG': 0.73813999999999991, 'CCT': 0.66512000000000004,
             'CTC': 0.66899999999999993, 'AGC': 0.81028000000000011, 'ACA': 0.66819999999999991,
             'AGA': 0.7396600000000001, 'CAT': 0.80928000000000011, 'AAT': 0.80779999999999996,
             'ATT': 0.73552000000000017, 'CTG': 0.47383999999999993, 'CTA': 0.47404000000000002,
             'ACT': 0.66971999999999998, 'CAC': 0.80871999999999999, 'ACG': 0.66917999999999989,
             'CAA': 0.80668000000000006, 'AGT': 0.80940000000000001, 'CAG': 0.80868000000000007,
             'CCG': 0.6644199999999999, 'CCC': 0.6670799999999999, 'TAT': 0.81141999999999992,
             'GGT': 0.66510000000000002, 'TGT': 0.80779999999999996, 'CGA': 0.59169999999999989,
             'CCA': 0.66690000000000016, 'CGC': 0.66578000000000004, 'GAT': 0.8096399999999998,
             'CGG': 0.59478000000000009, 'CTT': 0.66659999999999997, 'TGC': 0.80824000000000007,
             'GGG': 0.6641800000000001, 'TAG': 0.80645999999999995, 'GGA': 0.66880000000000006,
             'TAA': 0.61687999999999998, 'GGC': 0.66660000000000008, 'TAC': 0.80803999999999998,
             'GAG': 0.81069999999999998, 'TCG': 0.66408, 'TTA': 0.61682000000000003,
             'GAC': 0.80948000000000009, 'CGT': 0.6698400000000001, 'TTT': 0.80830000000000002,
             'TCA': 0.66878000000000004, 'GCA': 0.66555999999999993, 'GTA': 0.66721999999999992,
             'GCC': 0.6660600000000001, 'GTC': 0.66648000000000007, 'GCG': 0.66559999999999997,
             'GTG': 0.66546000000000005, 'TTC': 0.80668000000000006, 'GTT': 0.66598000000000002,
             'GCT': 0.66961999999999999, 'TGA': 0.80687999999999993, 'TTG': 0.61943999999999999,
             'TCC': 0.66835999999999995, 'TGG': 1.0, 'GAA': 0.80825999999999998,
             'TCT': 0.66670000000000007}

    _N_codons = len(sequence)/3

    #define recipient variables

    _dNdS = []
    _dNdS_par = []

    for _trial in range(N):#loop over number of trials

        _Ns = 0. #recipient for Nonsynonymous mutations

        #range of codon steps along the given sequence
        for _step in np.arange(0,len(sequence),3):
            #define the fished codon
            _fished_codon = sequence[_step:_step+3]
            #use the codon_mutate() to model mutation
            if method=='derived': _Ns+=codon_mutate_TiTv_GC(_fished_codon,
                                                            gc_content=gc,
                                                            ti=ti)
            #use a binomial random trial to determine if the mutation is Ns
            if method=='binomial': _Ns+=np.random.binomial(1, _pr_N[_fished_codon])

        _eNS = calculate_NS(sequence)
        _eN = 0.
        _eS = 0.
        for (_a,_b) in _eNS:
            _eN+=_a
            _eS+=_b
        try:
            _dNdS.append((_Ns/_eN)/((_N_codons-_Ns)/_eS))
            _dNdS_par.append((_Ns/_eN,(_N_codons-_Ns)/_eS))
        except:
            pass

    if pars: return _dNdS_par
    else: return _dNdS

def get_dNdS(annotated_dataset, gene_list=None, het_list=None, pars=False):
    """Calculate dNdS for a the annotated dataset.
    Same can be done to a subset of genes speicified
    by 'gene_list'.

    Assumes a specific data format - see below.

    ARGUMENTS:
    ----------
    annotated_dataset: annotation dict,
    example of value[0] in dict:
    ['Rv3885c', 'synonymous', '39', 'Val/V-Val/V', 'GTC-GTT']

    gene_list: list of genes of interest

    het_list: list of heterozygous positions of interest

    pars: return dN and dS separately - in case there are too
          mutations to effectively calculate dN/dS


    OUTPUT:
    -------
    dNdS: calculated dN/dS,float
    N: number of codons used for the calculation, float

    NOTES:
    ------
    First version. The format in which the data is inputed will
    probably change.

    Requires the calculate_NS() to be defined.
    """

    _concatenation = ''
    _N = 0.
    _S = 0.

    _dict_of_interest = {}

    if het_list!=None:
        for _het in het_list:
            _het = str(_het)
            if _het in annotated_dataset:
                _dict_of_interest[_het] = annotated_dataset[_het]

    if het_list==None:
        _dict_of_interest = annotated_dataset

    if gene_list==None:
        for _val in annotated_dataset.values():
            _info = _val[0]
            if _info[1]=='nonsynonymous':
                _N+=1.
                _codon = _info[4][:3]
                _concatenation+=_codon
            if _info[1]=='synonymous':
                _S+=1.
                _codon = _info[4][:3]
                _concatenation+=_codon

    if gene_list!=None:
        for _val in _dict_of_interest.values():
            _info = _val[0]
            if _info[1]=='nonsynonymous' and _info[0] in gene_list:
                _N+=1.
                _codon = _info[4][:3]
                _concatenation+=_codon
            if _info[1]=='synonymous' and _info[0] in gene_list:
                _S+=1.
                _codon = _info[4][:3]
                _concatenation+=_codon

    _NS = calculate_NS(_concatenation)

    _eN = 0.
    _eS = 0.

    for (_a,_b) in _NS:
        _eN+=_a
        _eS+=_b

    if pars: return _N/_eN, _S/_eS, _N+_S
    else: return ((_N/_eN)/(_S/_eS)), _N+_S

def pNS_CI(data, min_freq=0.005, N=1000, time=0, CI=95, codon_info='WT_CODON', parameters=False):
    """Calcualte pNS and their confidence intervals based on specific
    parameters.

    INPUTS:
    -------
    data: pd.DataFrame, like 'ALL'
    min_freq: float, minimum frequency for the allele
    N: int, number of bootstrapping iterations,
    time: 0|2|4|6|8, timepoint for which to calcualte
    CI: int, confidence interval
    codon_info: column label containing codon data
    parameters: bool, return pN and pS instead of pNS

    OUTPUTS:
    --------
    (pNS,CI_L, CI_H)

    NOTES:
    ------
    Depends on pandas as pd, scikits-learn utilities as sku,
    numpy as np, calculate_NS()
    """

    #define indices for empirical CI
    _CI_IND_L, _CI_IND_H = int((100-CI)/200.*N-1),int(-(100-CI)/200.*N)

    #define recipient variables
    _eNeS_results = []
    _oNoS_results = []

    #populate recipient variables based on parts above
    for trial in range(N):
        _data = sku.resample(data[(data['FREQUENCY']>=min_freq)&
                     (data['TIME']==time)])
        _eNeS_results.append(np.sum([[x[0],x[1]] for x in calculate_NS(''.join(_data[codon_info]))],0))
        _oNoS_results.append(np.sum(_data[['TYPE_NSY','TYPE_SYN']],axis=0))

    #reformat the output
    _eNeS_results = np.array(_eNeS_results,dtype=float)
    _oNoS_results = np.array(_oNoS_results,dtype=float)

    #calcualte the pN and pS as well as pNS
    _pN_pS = _oNoS_results/_eNeS_results
    _pNS = _pN_pS[:,0]/_pN_pS[:,1]

    #Get the specific indices for the measurements at CI edges
    _CI_INDS = np.argsort(_pNS)[[_CI_IND_L, _CI_IND_H]]

    if parameters:
        return np.median(_pN_pS, axis=0), _pN_pS[_CI_INDS[0]], _pN_pS[_CI_INDS[1]]

    else:
        return np.median(_pNS, axis=0), _pNS[_CI_INDS[0]], _pNS[_CI_INDS[1]]

#MARKOV CHAINS
def transition_matrix_PA(patient_matrix, annotation_dict, steps='Auto', min_freq=0, matrix_print=False, exclude=None,mutation_type=None):
    """Count the different types of transitions
    in our patient data and estimate transition
    probabilites.

    First convert the matrix to a binary matrix
    reflecting present/absent. The scoring depends
    on doubling the first column and subtracting
    the second. Possible outcomes are:

        0  - 0->0
        1  - 1->1
        2  - 1->0
        -1 - 0->1

    Incidentally the above is the order in which
    counts will emerge from "Counter".

    The transition probabilities are calculated as:

        p(0->0) = #0/(#0+#-1)
        p(1->1) = #1/(#1+#2)
        p(1->0) = #2/(#0+#-1)
        p(0->1) = #-1/(#1+#2)

    This calculation is done on the collated counts.

    INPUTS:
    ------
    patient_matrix: np.array <- populate_poly()
    annotation_dict: dict, {locus:['LOCUS_ID','synonymous/nonsynonymous'...]}
    steps: list, transitions to consider.
        steps[0] start transition, steps[1], number
        of steps to take.
        If 'Auto' it will take all the transitions
        for which there are data.
    min_freq: minimum frequency
    matrix_print: bool, print the transition matrix
    exclude: timepoints to exclude (e.g. missing data)
        Note that this index should reflect the
        postion within the patient_matrix.
    mutation_type ('synonymous'|'nonsynonymous'|'IGR'):
        only take mutations of a specific type

    OUTPUT:
    -------
    np.array(p(0->0), p(1->1), p(1->0), p(0->1))

    NOTES:
    ------

    This function is based on:
    http://www.stat.cmu.edu/~cshalizi/462/lectures/06/markov-mle.pdf
    """

    #Note that the patient matrix should be a
    #2D-np.array with locus at positon 1 and
    #then allele frequency in the following cols.

    #Start by transforming the array into a
    #Boolean array

    if exclude:
        patient_matrix = np.delete(patient_matrix,exclude,1)

    if steps=='Auto':
        #Just in case we want all the data for a patient
        #determine earliest timepoint with data
        _start = min(np.arange(1,8)[patient_matrix[:,1:].sum(0)>0])
        #determine the latest timepoint with data and calculate the # of steps
        n_steps = max(np.arange(1,8)[patient_matrix[:,1:].sum(0)>0])-_start
        #Turn steps into a list, where [0] is the start position and
        #[1] is the number of steps.
        steps = [_start,n_steps]

    _yesno = np.array(patient_matrix>min_freq,dtype=int)[:,np.arange(steps[0],1+np.sum(steps))]
    _result = _yesno[:,:-1]*2-_yesno[:,1:]

    if mutation_type: #filter by mutation type
       _ref = [float(k) for k,v in annotation_dict.items() if v[0][1]==mutation_type]
       _mtf = [_ind for _ind,_locus in enumerate(patient_matrix[:,0]) if _locus in _ref]
       _result = _yesno[_mtf,:-1]*2-_yesno[_mtf,1:]

    #Define count variables
    _c = np.array([Counter(_result.flatten()).get(x,0) for x in [0,1,2,-1]], dtype=float) #counts of 0->0, 1->1, 1->0, 0->1
    _t = np.array(_c+[_c[i] for i in [3,2,1,0]]) #transition counts to 0 and 1.

    if matrix_print:
        print('Transition matrix\n\t0\t1\n0\t%.3f\t%.3f\n1\t%.3f\t%.3f\n' %tuple((_c/_t)[i] for i in [0,2,3,1]))

    return _c/_t

def transition_matrix_PL(patient_matrix, annotation_dict, steps='Auto', min_freq=0, cutoff=0.01, matrix_print=False, exclude=None, mutation_type=None):
    """Count the different types of transitions
    in our patient data and estimate transition
    probabilites.

    First convert the matrix to a matrix
    reflecting absent(0)/present<=1%(1)/present>1%(2).
    The scoring depends on tripling the first column
    and subtracting the second.

    Possible outcomes are:

        0  - 0 -> 0
        1  - 1 -> 1+
        2  - 1 -> 1
        3  - 1 -> 0
        4  - 1+ -> 1+
        5  - 1+ -> 1
        6  - 1+ -> 0
        -1 - 0 -> 1
        -2 - 0 -> 1+

    Incidentally the above is the order in which
    counts will emerge from "Counter".

    The transition probabilities are calculated as:

        p(0->X) = #Xs/(#0+#-1+#-2)
        p(1->Y) = #Ys/(#1+#2+#3)
        p(1+->Z) = #Zs/(#4+#5+#6)

    This calculation is done on the collated counts.

    INPUTS:
    ------
    patient_matrix: np.array <- populate_poly()
    annotation_dict: dict, {locus:['LOCUS_ID','synonymous/nonsynonymous'...]}
    steps: list, transitions to consider.
        steps[0] start transition, steps[1], number
        of steps to take.
        If 'Auto' it will take all the transitions
        for which there are data.
    min_freq: minimum frequency
    cutoff: float, cutoff for inclusion into 1+
    matrix_print: bool, print the transition matrix
    exclude: timepoints to exclude (e.g. missing data)
        Note that this index should reflect the
        postion within the patient_matrix.
    mutation_type ('synonymous'|'nonsynonymous'|'IGR'):
        only take mutations of a specific type


    OUTPUT:
    -------
    np.array(p(0->0), p(1->1+), p(1->1), p(1->0),p(1+->1+), p(1+->1), p(1+->0), p(0->1+), p(0->1))

    NOTES:
    ------
    This script is based on:
    http://www.stat.cmu.edu/~cshalizi/462/lectures/06/markov-mle.pdf
    """

    #Note that the patient matrix should be a
    #2D-np.array with locus at positon 1 and
    #then allele frequency in the following cols.

    #Start by transforming the array into a
    #Boolean array

    if exclude:
        patient_matrix = np.delete(patient_matrix,exclude,1)

    if steps=='Auto':
        #Just in case we want all the data for a patient
        n_steps = list(np.array(patient_matrix>0,dtype=int).sum(axis=0)).index(0)-2
        steps = [1,n_steps] #Turn steps into a list, where [0] is the start position and [1] is the number of steps.

    _yesno = np.array(patient_matrix>min_freq,dtype=int)[:,np.arange(steps[0],1+np.sum(steps))]+np.array(patient_matrix>cutoff,dtype=int)[:,np.arange(steps[0],1+np.sum(steps))]
    _result = _yesno[:,:-1]*3-_yesno[:,1:]

    if mutation_type: #filter by mutation type
       _ref = [float(k) for k,v in annotation_dict.items() if v[0][1]==mutation_type]
       _mtf = [_ind for _ind,_locus in enumerate(patient_matrix[:,0]) if _locus in _ref]
       _result = _yesno[_mtf,:-1]*3-_yesno[_mtf,1:]

    #Define count variables
    _c = np.array([Counter(_result.flatten()).get(x,0) for x in [0,1,2,3,4,5,6,-2,-1]],dtype=float) #counts of 0->0, 1->1, 1->0, 0->1
    _t = np.array(_c+[_c[i] for i in [7,2,3,1,5,6,4,8,0]] + [_c[i] for i in [8,3,1,2,6,4,5,0,7]]) #transition counts to 0 and 1.

    if matrix_print:
        print('Transition matrix\n\t0\t1\t2\n0\t%.3f\t%.3f\t%.3f\n1\t%.3f\t%.3f\t%.3f\n2\t%.3f\t%.3f\t%.3f\n' %tuple((_c/_t)[i] for i in [0,3,6,-1,2,5,-2,1,4]))

    return _c/_t


def MC_PA_CI(patient_matrix, annotation_dict, min_freq=0, steps='Auto', exclude=None, mutation_type=None):
    """Merge transition data to re-sample and
    allow the calculation of transition
    probability CI.

    INPUTS:
    -------
    patient_matrix: np.array
    annotation_dict: dict, {locus:['LOCUS_ID','synonymous/nonsynonymous'...]}
    min_freq: minimum frequency threshold
    steps: list of steps to take into consideration
    exclude: define columns to exclude
    mutation_type ('synonymous'|'nonsynonymous'|'IGR'):
        only take mutations of a specific type

    OUTPUT:
    -------
    flattened np.array of simple outputs
    """

    if exclude: #exclude columns for which we don't have data
        patient_matrix = np.delete(patient_matrix,exclude,1)

    if steps=='Auto': #Just in case we want all the data for a patient
        #Just in case we want all the data for a patient
        #determine earliest timepoint with data
        _start = min(np.arange(1,8)[patient_matrix[:,1:].sum(0)>0])
        #determine the latest timepoint with data and calculate the # of steps
        n_steps = max(np.arange(1,8)[patient_matrix[:,1:].sum(0)>0])-_start
        #Turn steps into a list, where [0] is the start position and
        #[1] is the number of steps.
        steps = [_start,n_steps]

    _yesno = np.array(patient_matrix>min_freq,dtype=int)[:,np.arange(steps[0],1+np.sum(steps))]
    _result = _yesno[:,:-1]*2-_yesno[:,1:]

    if mutation_type: #filter by mutation type
       _ref = [float(k) for k,v in annotation_dict.items() if v[0][1]==mutation_type]
       _mtf = [_ind for _ind,_locus in enumerate(patient_matrix[:,0]) if _locus in _ref]
       _result = _yesno[_mtf,:-1]*2-_yesno[_mtf,1:]

    return _result.flatten()

def MC_PL_CI(patient_matrix, annotation_dict, steps='Auto', min_freq=0, cutoff=0.01, exclude=None, mutation_type=None):
    """Merge transition data to re-sample and
    allow the calculation of transition
    probability CI.

    INPUTS:
    -------
    patient_matrix: np.array
    annotation_dict: dict, {locus:['LOCUS_ID','synonymous/nonsynonymous'...]}
    steps: list of steps to take into consideration
    exclude: define columns to exclude
    mutation_type ('synonymous'|'nonsynonymous'|'IGR'):
        only take mutations of a specific type
    min_freq: minimum allele frequency

    OUTPUT:
    -------
    flattened np.array of simple outputs
    """

    if exclude:
        patient_matrix = np.delete(patient_matrix,exclude,1)

    if steps=='Auto':
        #Just in case we want all the data for a patient
        n_steps = list(np.array(patient_matrix>0,dtype=int).sum(axis=0)).index(0)-2
        steps = [1,n_steps] #Turn steps into a list, where [0] is the start position and [1] is the number of steps.

    _yesno = np.array(patient_matrix>min_freq,dtype=int)[:,np.arange(steps[0],1+np.sum(steps))]+np.array(patient_matrix>cutoff,dtype=int)[:,np.arange(steps[0],1+np.sum(steps))]
    _result = _yesno[:,:-1]*3-_yesno[:,1:]

    if mutation_type: #filter by mutation type
       _ref = [float(k) for k,v in annotation_dict.items() if v[0][1]==mutation_type]
       _mtf = [_ind for _ind,_locus in enumerate(patient_matrix[:,0]) if _locus in _ref]
       _result = _yesno[_mtf,:-1]*3-_yesno[_mtf,1:]

    return _result.flatten()

#EXCESS MUTATION ACCUMULATION
def selection_contingency_generator(genes,dataframe,exclude_check=False):
    """Generates a contingency table using a given set of
    genes. It will return two tables, one representing
    the status at time 0 and one covering the rest.

    INPUTS:
    -------
    genes: list, genes of interest
    dataframe: requires a df that has many of the
                same categories as ALL.
    exclude_check: False|list, check that all
                   problematic genes are removed.

    OUTPUTS:
    --------
    nd.array, nd.array
    the first array is the contingency table at t=0
    the second array is the contigency table at t!=0

    NOTES:
    ------
    Assumes the given data structure.

    Return array: [[NS_in, S_in],[NS_not, S_not]]

    EXAMPLE:
    --------
    ct0, ct1 = selection_contingency_generator(DR_set, ALL, exclude_check=excluded)

    """

    _df = dataframe.copy()

    if exclude_check:
        genes = [x for x in genes if x not in exclude_check]

    _df['INTEREST'] = [int(x in genes) for x in _df.GENE]

    _a = len(np.unique(_df.LOCUS[(_df.TYPE=='NSY')&(_df.INTEREST==1)&(_df.TIME==0)]))
    _a2 = len(np.unique(_df.LOCUS[(_df.TYPE=='NSY')&(_df.INTEREST==1)]))-_a
    _b = len(np.unique(_df.LOCUS[(_df.TYPE=='NSY')&(_df.INTEREST==0)&(_df.TIME==0)]))
    _b2 = len(np.unique(_df.LOCUS[(_df.TYPE=='NSY')&(_df.INTEREST==0)]))-_b
    _c = len(np.unique(_df.LOCUS[(_df.TYPE=='SYN')&(_df.INTEREST==1)&(_df.TIME==0)]))
    _c2 = len(np.unique(_df.LOCUS[(_df.TYPE=='SYN')&(_df.INTEREST==1)]))-_c
    _d = len(np.unique(_df.LOCUS[(_df.TYPE=='SYN')&(_df.INTEREST==0)&(_df.TIME==0)]))
    _d2 = len(np.unique(_df.LOCUS[(_df.TYPE=='SYN')&(_df.INTEREST==0)]))-_d

    return np.array([[_a,_c],[_b,_d]]), np.array([[_a2,_c2],[_b2,_d2]])

def geneset_size(gene_list, locus_dict):
    """Gives you the cumulative size (in bp) of the specified
    genes in the genome.

    INPUTS:
    -------
    gene_list: list, genes of interest
    locus_dict: dict, keys are genes in gene_list
                value is gene length float.

    OUTPUT:
    -------
    float, genomic size in bp.

    NOTES:
    ------
    This function only calculates the size of genes
    that are present in the locus_dict, if a gene is
    not present it will not contribute to the tally.
    """

    total_size = 0.

    for gene in gene_list:
        total_size+=locus_dict.get(gene,0)

    return total_size

def get_KEGG_pathway_genes(organism='mtu',pathway='mtu00010'):
    """Uses requests to query the KEGG REST API
    and retrieve all the genes in a given pathway.

    INPUTS:
    -------
    organism: str, KEGG organism code, mtu for H37Rv
    pathway: str, e.g. mtu00100 is H37Rv glycolysis

    OUTPUTS:
    --------
    list, genes in the pathway

    NOTES:
    ------
    Needs a connection and the requests module as rq
    """

    _query = 'http://rest.kegg.jp/link/%s/%s' %(organism, pathway)

    _fetch = rq.get(_query)
    _entries = _fetch.content.strip().split('\n')

    _genes = [x.split('\t')[-1].split(':')[-1] for x in _entries]

    return _genes

def excess_mutation(genes,dataframe,genome_size,locus_dict,exclude_check=False,H37Rv=True):
    """Performs multiple tests to determine
    whether or not the number of mutations in
    a locus is higher than expected.

    INPUTS:
    -------
    genes: list, genes of interest
    dataframe: requires a df that has many of the
                same categories as ALL.
    genome_size: float, denominator
    locus_dict: dict, keys are genes in gene_list
                value is gene length float.
    exclude_check: False|list, check that all
                   problematic genes are removed.
    H37Rv: boolean, whether or not to expect
                H37Rv identifiers.

    OUTPUTS:
    --------
    dict of data and statistics

    NOTES:
    ------
    This approach explicity ignores mutations that fall
    outside of coding regions for the sake of simplicity.
    Also, the counts are nonredundant (i.e. each allele
    is only counted once, even if it occurs at multiple
    timepoints.

    EXAMPLE:
    --------
    excess_mutation(DR_set, ALL, len(MTB_GENOME), SIZE_DICT, exclude_check=excluded)

    """

    _df = dataframe.copy()
    if H37Rv:
        genome_size = geneset_size([x for x in locus_dict.keys() if len(x) in [6,7]], locus_dict)

    if exclude_check:
        genes = [x for x in genes if x not in exclude_check]
        genome_size-=geneset_size(exclude_check, locus_dict)

    _p_hit = geneset_size(genes, locus_dict)/genome_size

    _df['INTEREST'] = [int(x in genes) for x in _df.GENE]

    _a = len(np.unique(_df.LOCUS[(_df.SNP_TYPE=='NSY')&(_df.INTEREST==1)&(_df.TIME==0)]))
    _a2 = len(np.unique(_df.LOCUS[(_df.SNP_TYPE=='NSY')&(_df.INTEREST==1)]))-_a
    _b = len(np.unique(_df.LOCUS[(_df.SNP_TYPE=='NSY')&(_df.INTEREST==0)&(_df.TIME==0)]))
    _b2 = len(np.unique(_df.LOCUS[(_df.SNP_TYPE=='NSY')&(_df.INTEREST==0)]))-_b
    _c = len(np.unique(_df.LOCUS[(_df.SNP_TYPE=='SYN')&(_df.INTEREST==1)&(_df.TIME==0)]))
    _c2 = len(np.unique(_df.LOCUS[(_df.SNP_TYPE=='SYN')&(_df.INTEREST==1)]))-_c
    _d = len(np.unique(_df.LOCUS[(_df.SNP_TYPE=='SYN')&(_df.INTEREST==0)&(_df.TIME==0)]))
    _d2 = len(np.unique(_df.LOCUS[(_df.SNP_TYPE=='SYN')&(_df.INTEREST==0)]))-_d

    _ct0, _ct1 = np.array([[_a,_c],[_b,_d]]), np.array([[_a2,_c2],[_b2,_d2]])

    #_N = len(np.unique(_df.LOCUS))
    _N = (_ct0+_ct1).sum()

    _result = {}

    #Enrichment of mutations following treatment initiation ie more mutations than expected once treatment has comenced
    try:
        _result['Enrichment'] = ss.chi2_contingency(np.vstack((_ct0.sum(1),_ct1.sum(1))),lambda_='log-likelihood')
    except:
        _result['Enrichment'] = ss.fisher_exact(np.vstack((_ct0.sum(1),_ct1.sum(1))))
    #Excess of mutations within given region (G-test, based on probability null)
    try:
        _result['Excess_Gtest'] = ss.chi2_contingency(np.vstack(((_ct0+_ct1).sum(1), np.array([_N*_p_hit, _N*(1-_p_hit)]))),lambda_='log-likelihood')
    except:
        _result['Excess_Fisher'] = ss.fisher_exact(np.vstack(((_ct0+_ct1).sum(1), np.array([_N*_p_hit, _N*(1-_p_hit)]))))
    #Excess of mutations within given region (Binomial tet, based on probability null)
    try:
        _result['Excess_binomial'] = ss.binom_test([_a+_a2+_c+_c2,_N-(_a+_a2+_c+_c2)],p=_p_hit,alternative='greater')
    except:
        _result['Excess_binomial'] = 'Failed'
    #Relative accumulation of nonsynonymous and synonymous polymorphisms
    try:
        _result['NS_binomial'] = ss.binom_test(_ct1[0],p=float(_ct1[1,0])/np.sum(_ct1[1]),alternative='greater')
    except:
        _result['NS_binomial'] = 'Failed'
    try:
        _result['Excess_binomial_treatment'] = ss.binom_test(sum(_ct1[0]), sum(_ct1),p=_p_hit,alternative='greater')
    except:
        _result['Excess_binomial_treatment'] = 'Failed'
    try:
        _result['NS_all'] = ss.chi2_contingency((_ct0+_ct1),lambda_='log-likelihood')
    except:
        _result['NS_all'] = ss.fisher_exact((_ct0+_ct1))
    #Relative accumulation of nonsynonymous and synonymous polymorphisms following treatment
    try:
        _result['NS_treatment'] = ss.chi2_contingency(_ct1,lambda_='log-likelihood')
    except:
        _result['NS_treatment'] = ss.fisher_exact(_ct1)
    #Relative accumulation of nonsynonymous and synonymous polymorphisms prior to treatment
    try:
        _result['NS_diagnosis'] = ss.chi2_contingency(_ct0,lambda_='log-likelihood')
    except:
        _result['NS_diagnosis'] = ss.fisher_exact(_ct0)
    #Data used to calculate statistics
    _result['Data'] = [_ct0, _ct1, _p_hit, _N, genome_size, geneset_size(genes, locus_dict)]

    return _result

def ScoreData(data,mapper,output='sum'):
    """Scores a collection of data based on a score mapper.

    INPUTS:
    -------
    data: 1darray-like
    mapper: dict, items in 'data' are keys
            values are int or float.
    output: 'sum'|'mean'

    OUTPUT:
    --------
    float

    NOTES:
    ------
    Expects the items in a list to be the same as keys in the
    mapper. Uses dict.get() returning 0 for queries that are
    not in the mapper dictionary.

    EXAMPLE:
    --------
    E.g. a list of codons ('codons') needs to be scored for
    the probability of being mutated to a nonsynonymous codon
    ('scores'):

    codons = ['ATG','GGT']
    scores = {'ATG': 1., 'GGT': 0.67}
    ScoreData(codons, scores, output='mean')

    >>> 0.835
    """
    _score = sum([mapper.get(_x,0) for _x in data])

    if output=='mean':
        return float(_score)/len(data)

    else:
        return float(_score)

def BinomialExcess(genes,dataframe,genome_size,locus_dict,score_dict,exclude_check=False,H37Rv=True,codons='WT_CODON'):
    """Performs binomal tests to determine
    whether or not the number of mutations in
    a locus is higher than expected and whether
    they are more likely than expected to be NSY.

    INPUTS:
    -------
    genes: list, genes of interest
    dataframe: requires a df that has many of the
                same categories as ALL.
    genome_size: float, denominator
    locus_dict: dict, keys are genes in gene_list
                value is gene length float.
    score_dict: dict, keys are codons in dataframe[codons]
                value is probabiliy of NSY float.
    exclude_check: False|list, check that all
                   problematic genes are removed.
    H37Rv: boolean, whether or not to expect
                H37Rv identifiers.
    codons: str, category in dataframe containing codon info.

    OUTPUTS:
    --------
    list(p_excess_mutations[2-tail], p_excess_NSY[2-tail],
         p_excess_mutations[1-tail], p_excess_NSY[1-tail])

    NOTES:
    ------
    This approach explicity ignores mutations that fall
    outside of coding regions for the sake of simplicity.
    Also, the counts are nonredundant (i.e. each allele
    is only counted once, even if it occurs at multiple
    timepoints. It also needs a dataframe where wild type
    codons are specified as it estimates the probability of
    observing a given number of NSY - diffrent models of
    substitution can be accounted for by using specific
    score_dict.

    EXAMPLE:
    --------
    BinomialExcess(DR_set, ALL, len(MTB_GENOME), SIZE_DICT, Pr_N, exclude_check=excluded)

    """
    _df = dataframe.copy()
    if H37Rv:
        genome_size = geneset_size([x for x in locus_dict.keys() if len(x) in [6,7]], locus_dict)

    if exclude_check:
        genes = [x for x in genes if x not in exclude_check]
        genome_size-=geneset_size(exclude_check, locus_dict)

    _p_hit = geneset_size(genes, locus_dict)/genome_size

    _df['INTEREST'] = [int(x in genes) for x in _df.GENE]

    _a2 = len(np.unique(_df.LOCUS[(_df.TYPE=='NSY')&(_df.INTEREST==1)&(_df.TIME!=0)]))
    _c2 = len(np.unique(_df.LOCUS[(_df.TYPE=='SYN')&(_df.INTEREST==1)&(_df.TIME!=0)]))
    _l2 = len(np.unique(_df.LOCUS[(_df['TYPE'].isin(['NSY','SYN']))&(_df.TIME!=0)]))

    _p_N = ScoreData(list(_df[_df['TYPE'].isin(['NSY','SYN'])].drop_duplicates('LOCUS')[codons]), score_dict, output='mean')

    return [[_a2,_c2,_l2],ss.binom_test((_a2+_c2),_l2,p=_p_hit),
            ss.binom_test(_a2,(_a2+_c2),p=_p_N),
            ss.binom.sf((_a2+_c2)-1, _l2, _p_hit),
            ss.binom.sf(_a2-1, (_a2+_c2), _p_N)]



#GENE ANNOTATION
def mutated_gene_id(SNP_position, gene_matrix):
    """Identify the location of a SNP.

    INPUTS:
    -------
    SNP_position: locus, integer
    gene_matrix: ND-array, start..stop

    OUTPUT:
    -------
    int, upstream edge of mutated locus.

    NOTES:
    ------
    This was established to be used on our MTB_anc
    annotation. It will return the start positon of
    a gene to be used with REF_annotation.

    REF_annotation={start:[end, orientation, locus_ID],...}

    The rationale is: the only case where (start-SNP)*(end-SNP)
    is negative will be when start<SNP<end. Both start<end<SNP
    and SNP<start<end will produce positive values when multiplied.

    NB. If the SNP falls exactly at the edge of a locus the product
    will be zero. Therefore we're looking for non-positive products.
    """

    return int(gene_matrix[np.product(gene_matrix-SNP_position, axis=1)<=0][0][0])

def codon_fetch(SNP_position, gene_cooridnates, genome, return_params=False):
    """Extracts the codon based on a SNP coordinate.

    INPUTS:
    -------
    SNP_position: locus, integer
    gene_coordinate: tuple, pythonic indices
    genome: string, fasta of the genome
    return_params: bool, return position information

    OUTPUT:
    -------
    codon: if coding, else IGR

    NOTES:
    ------
    This was established to be used on our MTB_anc
    annotation
    """

    _gene_start = gene_cooridnates[0]
    _gene_end = gene_cooridnates[1]

    if _gene_start<_gene_end: #Forward gene
        _into_gene = SNP_position-_gene_start
        _direction = 1 #set direction, 1 for F, -1 for R
    if _gene_start>_gene_end: #Reverse gene
        _into_gene = _gene_start-SNP_position
        _direction = -1

    _codon_position = int(_into_gene)/3+1 #codon
    _position_in_codon = (int(_into_gene)%3+1)*_direction #base in codon [1,2,3] for F or [-1,-2,-3] for R

    _slice = {1:(0,3), 2:(-1,2), 3:(-2,1), -1:(-2,1), -2:(-1,2), -3:(0,3)} #define codon locations

    _begin = SNP_position+_slice[_position_in_codon][0] #where to begin slicing a genome
    _end = SNP_position+_slice[_position_in_codon][1] #where to end slicing a genome

    if _direction==1:
        if return_params: return genome[_begin:_end], _codon_position, _position_in_codon
        else: return genome[_begin:_end]

    if _direction==-1:
        #Need the reverse complement of the codon. Used the translate function from
        #the string module.
        xmz = 'ATGC' #the basis for the replacement
        mzx = 'TACG'
        complement=(genome[_begin:_end].translate(str.maketrans(xmz, mzx))) #generate complement

        if return_params: return complement[::-1], _codon_position, _position_in_codon
        else: return complement[::-1] #reverse complement

def mutation_classifier(reference_base, alternative_base):

    """
    Determine mutation type: transition vs transversion

    INPUTS:
    -------
    refrence_base: A|C|G|T reference base, upper case string
    alternative_base: A|C|G|T mutated base, upper case string

    OUTPUT:
    -------
    int: 0 for transversion, 1 for transition.

    NOTES:
    ------
    see: https://en.wikipedia.org/wiki/Transition_(genetics)
    """

    _mutation_class = {'A':'G', 'C':'T', 'G':'A', 'T':'C'}
    return int(_mutation_class[reference_base]==alternative_base)

def SNP_annotate(SNP_position, wt_base, mut_base, locus_matrix, reference, genome, codon_table):

    """Generate a basic annotation for SNPs.

    INPUTS:
    -------
    SNP_position: locus, integer
    wt_base: str, reference base
    mut_base: str, mutant base
    locus_matrix: 2D-array, [[start,end],...]
    reference: dict, {start:[end,direction,gene],...}
    genome: genome sequence in FASTA
    codon_table: {codon: AA,...}

    OUTPUT:
    -------
    list: [Rv, mutation_type, codon positon, AA/AA, W_codon-M_codon, Ti]
        if IGR: [Rv-Rv, 'IGR', -, ---, _N-M_, Ti]
        if RNA: [name, 'RNA', '-', '---', W#posM, Ti]

    NOTES:
    ------
    This was established to be used on our MTB_anc
    annotation and with the H37Rv_QC snp input data:
        pos wt_base alt_base    freq
    """

    _gene_internal_id = mutated_gene_id(SNP_position, locus_matrix)
    _gene = reference[_gene_internal_id][2]

    if reference[_gene_internal_id][1]=='+':
        _start, _end = _gene_internal_id,reference[_gene_internal_id][0] #define start and end of CDS
        _wt_codon, _codon_pos, _pos_in_codon = codon_fetch(SNP_position, (_start,_end), genome, return_params=True) #get codon information
        _mut_codon = _wt_codon[:abs(_pos_in_codon)-1]+mut_base+_wt_codon[abs(_pos_in_codon):] #derive the mutated codon, change the position affected by the mutation keeping the other two
        #Call the mutation type
        if codon_table[_wt_codon]==codon_table[_mut_codon]:
            _mutation_type='synonymous'
        if codon_table[_wt_codon]!=codon_table[_mut_codon]:
            _mutation_type='nonsynonymous'
        #Generate the information output
        return [_gene, _mutation_type, _codon_pos, '%s/%s' %(codon_table[_wt_codon], codon_table[_mut_codon]), '%s-%s' %(_wt_codon, _mut_codon),mutation_classifier(wt_base, mut_base)]

    if reference[_gene_internal_id][1]=='-':
        _complementary_base = {'A':'T', 'C':'G', 'G':'C', 'T':'A'} #Assumes all SNP calls are done in the + orientation.
        _start, _end = reference[_gene_internal_id][0]-1, _gene_internal_id+2 #define start and end of CDS, need to tweak a bit, to make sure we stay in frame. Keep an eye out.
        #The values used here were determined based on a manual insepction of one gene.
        _wt_codon, _codon_pos, _pos_in_codon = codon_fetch(SNP_position, (_start,_end), genome, return_params=True) #get codon information
        _mut_codon = _wt_codon[:abs(_pos_in_codon)-1]+_complementary_base[mut_base]+_wt_codon[abs(_pos_in_codon):] #derive the mutated codon, change the position affected by the mutation keeping the other two
        #Call the mutation type
        if codon_table[_wt_codon]==codon_table[_mut_codon]:
            _mutation_type='synonymous'
        if codon_table[_wt_codon]!=codon_table[_mut_codon]:
            _mutation_type='nonsynonymous'
        #Generate the information output
        return [_gene, _mutation_type, _codon_pos, '%s/%s' %(codon_table[_wt_codon], codon_table[_mut_codon]), '%s-%s' %(_wt_codon, _mut_codon),mutation_classifier(wt_base, mut_base)]

    if reference[_gene_internal_id][1]=='I':
        _gene = reference[_gene_internal_id][2].split('_')[1] #Modify _gene to reflect only the position of the intergenic region.
        _len_to_symbol = {3:'+',6:'+', 7:'-', 9:'..'} #define dict to be used for the construction of the IGR tag string. NB this trick does not work for CCDC5079. In that case '..' is added on both sides of the gene...
        _up,_down = _len_to_symbol.get(len(_gene.split('-')[0]),'..'), _len_to_symbol.get(len(_gene.split('-')[1]),'..') #get +/- for F and R genes respectively
        _IGR_tag = '%s%d-%d%s' %(_up,SNP_position-_gene_internal_id,reference[_gene_internal_id][0]-SNP_position,_down) #generate IGR tag: +/- denotes whether the upstream or downstream if F or R. Distance from start/end of neighbouring genes
        return [_gene, 'IGR', '-', '---', _IGR_tag, mutation_classifier(wt_base, mut_base)] #define output list

    if reference[_gene_internal_id][1]=='R':
        _start, _end = _gene_internal_id,reference[_gene_internal_id][0] #get start and stop for ncRNA
        _into_gene = SNP_position-_start #get position in ncRNA
        _RNA_tag = '%s%d%s' %(wt_base, _into_gene, mut_base)
        return [_gene, 'RNA', '-', '---', _RNA_tag, mutation_classifier(wt_base, mut_base)] #define output list

#Define functions
def logistic(x, K, r, N0):
    """Generate a logistic curve

    Input
    -------
    x: 1D array
    K: carrying capacity
    r: Growth rate
    N0: starting population

    Output
    -------
    Calculates the population at time x given the parameters using a simple logistic model
    adapted from wikipedia: http://en.wikipedia.org/wiki/Logistic_function

    Notes
    -------
    Made to be used with scipy.optimize.curve_fit
    """
    return (K*N0*np.exp(r*x))/(K+N0*(np.exp(r*x)-1))

def exp_decay(x, N, k):
    """Calculate the exponential decay of a number.

    Input
    -------
    x: 1D array
    N: Starting population
    k: exponential decay constant

    Output
    -------
    1D array

    Notes
    -------
    Adapted form: http://en.wikipedia.org/wiki/Exponential_decay
    """
    return N*np.exp(-x/k)

def demo_dynamic(g, N0, r, K, k, point_A, point_S, r_=None, K_=None):
    """Demographic function of growth, cell death, sampling and
    regrowth.

    INPUT:
    ------
    g: number of generations to be conisdered
    N0: starting population (logistic)
    r: growth rate (logistic)
    K: carrying capacity (logistic)
    k: decay rate (exponential decay)
    point_A: generation of effective treatment
    point_S: generation of sampling
    r_: growth rate of re-growth
    K_: carrying capacity of re-growth

    OUTPUT:
    -------
    Population size dynamics

    NOTES:
    ------
    Modelled on:
    http://simupop.sourceforge.net/manual_release/build/userGuide_ch8_sec2.html
    therefore expected to be called as:
        demo_func = demo_dynamic(101, 100., np.log(2), 100000., 15., 30, 60)
        demo_func(70) : population size at generation 70

    """

    if r_==None: r_ = r
    if K_==None: K_ = K

    _unperturbed = logistic(np.arange(g), K, r, N0)
    _decayed = exp_decay(np.arange(g), _unperturbed[point_A], k)
    _regrown = logistic(np.arange(g), K_, r_, _decayed[point_S-point_A]*0.001) #1/1000 dilution

    _overall = np.concatenate((_unperturbed[:point_A], _decayed[:point_S-point_A], _regrown))

    def snipett(gen):
        return _overall[gen]

    return snipett
