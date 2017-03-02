# Processed data
These data are provided in the data repository:
doi:10.5281/zenodo.322377

The files include:
- ALLELE_data.csv
- ALLELE_info.pkl
- CLINICAL_info.csv 	
- PHYLOGENY_strains.nwk

##Patient data
All clinical patient data are given in "PATIENT_data.csv" on the repository.
The categories in this file are coded as follows:
- *DRUG*: Drug susceptibility profile of patient
- *PATIENT_ID*: Patient ID
- *TIME*: time post treatment commencement (weeks)
- *TIME_TO_POSITIVITY_MGIT*: Week at which MTBC populations were first harvested from MGIT tubes
- *SEQUENCED_CULTURE*: Source of sequenced population (Culture conditions)
- *DEPTH*: Average read depth for the sample
- *NON_EFFICACIOUS*: Boolean classifier reflecting treatment efficacy (4+ drugs: 0, <4 drugs: 1)
- *vSNP_COUNT*:  number of v-SNPs at timepoint
- *fSNP_COUNT*:  number of f-SNPs at timepoint
- *EXPECTED_NSY*:  expected number of nonsynonymous mutations based on mutated codons
- *EXPECTED_SYN*:  expected number of synonymous mutations based on mutated codons
- *NEUTRAL_NSY*:  simulated number of nonsynonymous mutations based on given mutated codons under the assumption of neutrality
- *NEUTRAL_SYN*:  simulated number of synonymous mutations based on given mutated codons under the assumption of neutrality
- *OBSERVED_NSY*:  observed number of nonsynonymous mutations
- *OBSERVED_SYN*:  observed number of synonymous mutations
- *pNS*:  calculated pNS from observed mutations (empty if insufficient data)
- *pNS_NEUTRAL*:  calculated pNS from mutations simulated by assumption of neutrality

##v-SNP data
All collected and calculated data for all v-SNPs are given in "ALLELE_data.csv" on the repository.
The categories in this file are coded as follows:
- *PATIENT_ID*: Patient ID
- *RESISTANCE*:  Drug susceptibility profile of patient
- *FREQUENCY*: v-SNP frequency
- *LOG_FREQ*: log2 of v-SNP frequency
- *GENE*: ID tag of gene affected by v-SNP (H37Rv, CDS only)
- *LOCUS*:  genomic locus of v-SNPs
- *TRANSITION*: Boolean classifier reflecting weather the mutation is a transition
- *WT_CODON*: wild type codon affected by mutation
- *TIME*: time post treatment commencement (weeks)
- *SNP_TYPE*: v-SNP type
- *NON_EFFICACIOUS*: Boolean classifier reflecting treatment efficacy (4+ drugs: 0, <4 drugs: 1)
- *RECURRENT*: Boolean classifier reflecting whether a v-SNP is recurrent

##f-SNP data
All collected and calculated data for all f-SNPs are given in "FIXED_ALLELE_data.csv" on the repository.
The categories in this file are coded as follows:
- *PATIENT_ID*: Patient ID
- *RESISTANCE*: Drug susceptibility profile of patient
- *GENE*: ID tag of gene affected by v-SNP (CCDC5079, CDS only)
- *LOCUS*:  genomic locus of v-SNPs
- *TRANSITION*: Boolean classifier reflecting weather the mutation is a transition
- *WT_CODON*: wild type codon affected by mutation
- *SNP_TYPE*: v-SNP type
- *NON_EFFICACIOUS*: Boolean classifier reflecting treatment efficacy (efficacious: 0, non-efficacious: 1)
