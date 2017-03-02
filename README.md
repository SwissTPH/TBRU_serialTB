# TBRU_serialTB
Population dynamics of _Mycobacterium tuberculosis_ populations within hosts during treatment.


[![DOI](https://zenodo.org/badge/75634041.svg)](https://zenodo.org/badge/latestdoi/75634041)

The data pertinent to the analyses can be on zenodo, with DOI-references:
[![DOI](https://www.zenodo.org/badge/DOI/10.5281/zenodo.322377.svg)](https://doi.org/10.5281/zenodo.322377)

### The study was performed by:
Andrej Trauner\*, Qingyun Liu\*, Laura E. Via, Xin Liu, Xianglin Ruan, Lili Liang, Huimin Shi, Ying Chen, Ziling Wang, Ruixia Liang, Wei Zhang, Wang Wei, Jingcai Gao，Gang Sun, Daniela Brites, Kathleen England, Guolong Zhang, Sebastien Gagneux, Clifton E. Barry III, Qian Gao

\* equal contribution.

## Abstract for The within-host population dynamics of _Mycobacterium tuberculosis_ vary with treatment efficacy.

### Background:
Combination therapy is one of the most effective tools for limiting the emergence of drug resistance. Despite the widespread adoption of combination therapy across diseases, drug resistance rates continue to rise, leading to failing treatment regimens. The mechanisms underlying treatment failure are well studied, but the processes governing successful combination therapy are poorly understood. We addressed this question by studying the population dynamics of _Mycobacterium tuberculosis_ within tuberculosis patients undergoing treatment with different combinations of antibiotics.
### Results:
By combining very deep whole genome sequencing (~1,000-fold genome-wide coverage) with sequential sputum sampling, we were able to detect transient genetic diversity driven by the apparently continuous turnover of minor alleles, which could serve as the source of drug-resistant bacteria. However, we report that treatment efficacy had a clear impact on the population dynamics: sufficient drug pressure bore a clear signature of purifying selection leading to apparent genetic stability. In contrast, _M. tuberculosis_ populations subject to less drug pressure showed markedly different dynamics, including cases of acquisition of additional drug resistance.
### Conclusions:
Our findings show that for a pathogen like _M. tuberculosis_, which is well adapted to the human host, purifying selection constrains the evolutionary trajectory to resistance in effectively treated individuals. Nonetheless, we also report a continuous turnover of minor variants, which could give rise to the emergence of drug resistance in cases of drug pressure weakening. Monitoring bacterial population dynamics could therefore provide an informative metric for assessing the efficacy of novel drug combinations.


## Project structure
```
├── LICENSE
├── README.md
├── data #This is just a place-holder. All data available from Zenodo
│   ├── README.md #See this file for a better description
│   ├── external
│   ├── interim
│   ├── processed
│   └── raw
├── notebooks #jupyter notebooks for figures and analysis
├── reports #manuscript and figures
│   └── figures
├── requirements.txt
└── scripts #analysis scripts
    ├── variant_analysis #Downstram analysis
    │   ├── 1_Markov_chain.py
    │   ├── 2_pNS.py
    │   ├── 3_enrichment.py
    │   └── 4_demographic_simulations.py
    ├── variant_extraction #Variant calling scripts.
    │   ├── 0_Serial_DeepSeq_Variant_calling.sh
    │   ├── support_data
    │   │   ├── 1_Tuberculist.info
    │   │   ├── 2_genetic_codes
    │   │   ├── 3_H37Rv.fasta
    │   │   ├── 4_H37Rv_PPE_INSERT_loci_list
    │   │   └── PE_PGRS_list
    │   └── support_scripts
    │       ├── 0_sanger_triming_paired_ended.sh
    │       ├── 1_Mapping_bwa_paired_ended.sh
    │       ├── 1_format_trans.pl
    │       ├── 3.1_mix_pileup_merge.pl
    │       ├── 3_mix_extract.pl
    │       ├── 3a_mix_extract.pl
    │       ├── 5_lable_filter.py
    │       ├── 6_H37Rv_annotate.pl
    │       ├── tp_fp_test_new_151_201503.R
    │       └── varscan_lofreq_compare.pl
    └── variant_processing #Variant call processing scripts (python)
        ├── 1_parse_patient_FINAL.py
        ├── 2_process_patient_data.py
        ├── 3_allele_dataframe.py
        ├── 4_allele_data_polishing.py
        ├── 4_patient_dataframes.py
        ├── parse_ART_simulations.py
        ├── parse_single_colony.py
        └── serial_functions.py
```
