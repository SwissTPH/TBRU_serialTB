#!/bin/bash

################################################################################
#
# Script was written for SNP calling from deep sequencing reads.
#
# (Primitive) usage: 
#	
#	./0_Serial_DeepSeq_Variant_calling.sh path/to/my.sorted.bam
#
# The script is expecting sorted BAM files, derived from trimmed fastq files.
#
#
# Outputs filtered SNP lists.
#
# Qingyun Liu & Andrej Trauner, 01. December 2016
#
################################################################################
#	NOTE
################################################################################
#
# We used the following approach to generate sorted BAM files.
#
#-------------------------------------------------------------------------------
# Trim paired-end fastq files
# ./support_scripts/0_sanger_triming_paired_ended.sh
#		Enter the raw FASTAQ file Directory: #give the path of fastq files
#		Are they using Sanger scoring?       #give the quality format (yes/no)
#
#-------------------------------------------------------------------------------
# mapping to template with bwa
# ./support_scripts/1_Mapping_bwa_paired_ended.sh
#		enter the .fq sequence Directory:  #give the path of trimmed fastq files
#		enter the Template Directory:      #give the path where genome template locates
#
#-------------------------------------------------------------------------------
# sort the sequences in bam file
# samtools sort my.bam -o my.sort.bam
#
################################################################################


########################### Requirements ########################################
#perl/5.18.2
#python/2.7.10
#R/3.2.1
#SAMtools/1.2
#VarScan/2.3.6
#LoFreq/0.6.1 (java)

########################### Definitions ########################################
bam=$1

FILENAME=$(basename $1) ### take everything from the script argument that comes after the slash 
BASENAME=$(echo ${FILENAME%.sort.bam}) ### take everything from the script argument that comes after the slash and add '.fastq'
PATHNAME=$(dirname $1) ### take everything from the script argument that comes before the slash 
REFERENCE=./support_data/3_H37Rv.fasta

VARSCAN= #/path/to/VarScan.v2.3.6.jar
LOFREQ_SNP= #/path/to/lofreq_snpcaller.py
LOFREQ_FILTER= #/path/to/lofreq_filter.py

echo "processing $PATHNAME/$FILENAME"
echo "Generating mpileup for $BASENAME"

#This step generates a mpileup and filters out base_q<30, mapping_q<20, and outputs position data
#the pileup file is stored in temporary storage to be deleted at the end of the run.

samtools mpileup -q 30 -Q 20 -BOf $REFERENCE $1 > $BASENAME.pileup

echo "mpileup for $BASENAME made"

#Use LoFreq to filter out SNPs that:
# don't comply with the default quality parameters of LoFreq.
# have a coverage of less than 50, or greater than 3000

echo "Calling variants in $BASENAME with LoFreq"
#call low frequency SNPs using lofreq
python $LOFREQ_SNP -f $REFERENCE -b $BASENAME.sort.bam -o $BASENAME.lofreq.snp

#filter false positive SNPs with lofreq
python $LOFREQ_FILTER --strandbias-holmbonf --min-cov 50 --max-cov 3000 -i $BASENAME.lofreq.snp -o $BASENAME.lofreq.filt.snp
echo "LoFreq finished"

#Use VarScan to filter out SNPs that have:
# fewer than 4 reads supporting them,
# a coverage of less than 50,
# min average quality <20,
# min frequency <0.5%
# max frequency <90%

echo "Calling variants in $BASENAME with VARSCAN"
java -jar $VARSCAN mpileup2snp $BASENAME.pileup --min-coverage 50 --min-reads2 4 --min-avg-qual 20 --min-var-freq 0.005 --min-freq-for-hom 0.9 --p-value 99e-02 --strand-filter 1 > $BASENAME.varscan.var
echo "VARSCAN finished"

#Merge the LoFreq and Varscan output

echo "Merge congruent variants in $BASENAME"
perl /support_scripts/varscan_lofreq_compare.pl $BASENAME.lofreq.filt.snp $BASENAME.varscan.var > $BASENAME.var
echo "Merge done"

#Re-format the file for downstream processing

echo "Re-formatting the file"
perl /support_scripts/1_format_trans.pl $BASENAME.var > $BASENAME.for
echo "Re-formatting done"

#Get heterozygous SNP candidates:

echo "Extracting low frequency SNP candidates"
perl /support_scripts/3_mix_extract.pl $BASENAME.for > $BASENAME.lofreq
echo "Extraction done"

#Merge info to produce the data to use for a Kolmogorov-Smirnov test of distributions:

echo "Merging the pileups and low freq annotations"
perl /support_scripts/3.1_mix_pileup_merge.pl $BASENAME.pileup $BASENAME.lofreq > $BASENAME.merge
echo "Merging done"

#R script to remove SNPs whose support comes from the tail-end of the distribution.

echo "Running R: Kolmogorov–Smirnov test of SNP support"
Rscript /support_scripts/tp_fp_test_new_151_201503.R $BASENAME.merge
echo "KS test done"

#Python script to remove SNPs at specific sites, if frequency is less than 0.5% and the KS
#test is positive (p<0.05).

echo "Running Python: filter out SNPs from problematic regions, KS+ SNPs"
python /support_scripts/5_lable_filter.py $BASENAME.merge.lable /support_data/4_H37Rv_PPE_INSERT_loci_list
echo "$BASENAME sample fully processed"

echo "Annotating filtered SNP calls"
#annotate SNPs using H37Rv template(NC_000962.2)
perl H37Rv_annotation.pl /support_data/1_Tuberculist.info /support_data/2_genetic_codes $REFERENCE $BASENAME.filter.snp > $BASENAME.ano.snp
echo "Annotation done"