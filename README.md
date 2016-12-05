# TBRU_serialTB
Population dynamics of Mycobacterium tuberculosis within hosts during treatment.

Repo to be updated soon.

## Project structure
```
├── LICENSE
├── README.md
├── data
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
    │   └── 3_enrichment.py
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
