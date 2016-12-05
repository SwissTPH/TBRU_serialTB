#!/bin/bash
#!/usr/bin/perl

#before mapping the sequqnces with template, please make sure that the sequence quality format is in Sanger style
#And the .fq files have been trimmed
echo -n "enter the .fq sequence Directory: "
read fqd
fqd=${fqd%/}/

echo -n "enter the Template Directory: "
read tmpd
tmpd=${tmpd%/}/


#acquire paired .fq files for paired end mapping
m=0
n=0
l=0
for fq in `ls $fqd*.fq`
do
  fq_name_ori[$m]=${fq#$fqd}
  fq_a[$m]=${fq_name_ori[$m]#*_*_}
  fq_b[$m]=${fq_a[$m]%_*_*}

  if [ "${fq_b[$m]}" = "1" ];
  then
    fq1[$n]=$fq
    n=`expr $n + 1`
  fi

  if [ "${fq_b[$m]}" = "2" ];
  then
    fq2[$l]=$fq
    l=`expr $l + 1`
  fi
m=`expr $m + 1`
done

j=0
k=0
while [ $j -lt $n ]
do
  h=0
    while [ $h -lt $l ]
    do
      if [ "${fq1[$j]%_?_*_*}" = "${fq2[$h]%_?_*_*}" ];
      then
        fq1_name[$k]=${fq1[$j]#$fqd}
        fq2_name[$k]=${fq2[$h]#$fqd}
        fq1_path[$k]=${fq1[$j]}
        fq2_path[$k]=${fq2[$h]}
        fq_single_path[$k]=${fq1[$j]%_?_*_*}_single_san_trimmed.fq
        fq_name[$k]=${fq1_name[$j]%_?_*_*}
        k=`expr $k + 1`
        break
      else
        h=`expr $h + 1`
      fi
    done
j=`expr $j + 1`
done

echo -e "\n $k paired *.fq file(s) found! \n"


#acquire template file
  b=0
  for tmp in `ls $tmpd*.bwt | sed -e "s:$tmpd::;s/.bwt$//"`
  do
  tmpfile[$b]=$tmp
  echo " `expr $b + 1`. ${tmpfile[$b]}"
  b=`expr $b + 1`
  done

if [ $b -gt 1 ];then
  echo " default: 1. ${tmpfile[0]}"
  echo -n " Multiple templates found! Please select one: "
  read tmpnum

  if [ "$tmpnum" != "" ];then
    if [ $tmpnum -le $b ] && [ $tmpnum -gt 0 ];then
      template=${tmpfile[`expr $tmpnum - 1`]}
      echo -e "\n Template file is $template"
    fi
  else
      template=${tmpfile[0]}
      echo -e "\n Template file is $template"
  fi

else
  template=${tmpfile[0]}
  echo -e "\n Template file is $template"
fi



#check the .fai file
for fai in `ls template/fai/*.fai | sed -e "s:template/fai/::"`
do
  if [ "${fai%.fai}" = "$template" ];
  then
    faifile=$fai
    break
  else
    faifile="none"
  fi
done

if [ "$faifile" = "none" ];
then
  echo ".fai file for $template not found! Please creat one!"
exit
fi


#perform bwa paired end mapping
echo -e "\nBwa mapping started ...\n"
c=0
while [ $c -lt $k ]
do
echo -e "1. Initial mapping of strain ${fq_name[$c]} by bwa ...\n"
saifile1[$c]=${fq_name[$c]}.1.$template.bwa.sai
saifile2[$c]=${fq_name[$c]}.2.$template.bwa.sai
saifile_single[$c]=${fq_name[$c]}.single.$template.bwa.sai

paired_samfile[$c]=${fq_name[$c]}.$template.bwa.paired.sam
single_samfile[$c]=${fq_name[$c]}.$template.bwa.single.sam
samfile[$c]=${fq_name[$c]}.$template.bwa.sam

#test######
#echo -e "fq1_path= ${fq1_path[$c]}\n"
#echo -e "fq2_path= ${fq2_path[$c]}\n"
#echo -e "sam_ori= ./output/sam/${paired_samfile[$c]}\n"
###########

bwa aln -R 1 $tmpd$template ${fq1_path[$c]} > ./output/sam/${saifile1[$c]}
bwa aln -R 1 $tmpd$template ${fq2_path[$c]} > ./output/sam/${saifile2[$c]}
bwa aln -R 1 $tmpd$template ${fq_single_path[$c]} > ./output/sam/${saifile_single[$c]}

bwa sampe -a 1000 -n 1 -N 1 $tmpd$template ./output/sam/${saifile1[$c]} ./output/sam/${saifile2[$c]} ${fq1_path[$c]} ${fq2_path[$c]} > ./output/sam/${paired_samfile[$c]}
bwa samse -n 1 $tmpd$template ./output/sam/${saifile_single[$c]} ${fq_single_path[$c]} > ./output/sam/${single_samfile[$c]}


echo -e "\n2. SAM processing of strain ${fq_name[$c]%_?_*}\n"

echo -e "a) SAM -> BAM\n"
#SAM -> BAM
samtools view -bhSt template/fai/$faifile ./output/sam/${paired_samfile[$c]} -o samtools/bam/${paired_samfile[$c]%.sam}.bam
samtools view -bhSt template/fai/$faifile ./output/sam/${single_samfile[$c]} -o samtools/bam/${single_samfile[$c]%.sam}.bam
samtools merge samtools/bam/${samfile[$c]%.sam}.bam samtools/bam/${paired_samfile[$c]%.sam}.bam samtools/bam/${single_samfile[$c]%.sam}.bam
samtools sort samtools/bam/${samfile[$c]%.sam}.bam samtools/sort.bam/${samfile[$c]%.sam}.sort
samtools index samtools/sort.bam/${samfile[$c]%.sam}.sort.bam

echo -e "b) Mark duplicates\n"
#Mark duplicates
java -Xmx2g -jar /opt/bio/picard-tools/MarkDuplicates.jar INPUT=samtools/sort.bam/${samfile[$c]%.sam}.sort.bam OUTPUT=samtools/sort.bam/${samfile[$c]%.sam}.sort.rmdup.temp.bam REMOVE_DUPLICATES=true ASSUME_SORTED=true OPTICAL_DUPLICATE_PIXEL_DISTANCE=10 METRICS_FILE=samtools/sort.bam/${samfile[$c]%.sam}.sort.rmdup.metrix VALIDATION_STRINGENCY=SILENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000
#samtools index samtools/sort.bam/${samfile[$c]%.sam}.sort.rmdup.temp.bam

java -Xmx2g -jar /opt/bio/picard-tools/AddOrReplaceReadGroups.jar INPUT=samtools/sort.bam/${samfile[$c]%.sam}.sort.rmdup.temp.bam OUTPUT=samtools/sort.bam/${samfile[$c]%.sam}.sort.rmdup.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT RGLB=8 RGPL=illumina RGPU=1 RGSM=${samfile[$c]%.sam}
samtools index samtools/sort.bam/${samfile[$c]%.sam}.sort.rmdup.bam

echo -e "c) Indel realignment\n"
#Indel realignment
java -Xmx2g -jar /opt/bio/GATK/dist/GenomeAnalysisTK.jar -I samtools/sort.bam/${samfile[$c]%.sam}.sort.rmdup.bam -R ./template/fasta/${faifile%.fai} -T RealignerTargetCreator -maxInterval 20000 -o ./output/recalibration/${samfile[$c]%.sam}.sort.rmdup.indel.intervals
java -Xmx4g -jar /opt/bio/GATK/dist/GenomeAnalysisTK.jar -I samtools/sort.bam/${samfile[$c]%.sam}.sort.rmdup.bam -R ./template/fasta/${faifile%.fai} -T IndelRealigner -targetIntervals ./output/recalibration/${samfile[$c]%.sam}.sort.rmdup.indel.intervals -o samtools/sort.bam/${samfile[$c]%.sam}.sort.rmdup.realn.bam
#samtools index samtools/sort.bam/${samfile[$c]%.sam}.sort.rmdup.realn.bam


<<mark
echo -e "d) Base quality recalibration\n"
#Base quality recalibration
samtools mpileup -uf ./template/fai/${faifile%.fai} samtools/sort.bam/${samfile[$c]%.sam}.sort.rmdup.realn.bam > samtools/bcf/${samfile[$c]%.sam}.sort.rmdup.realn.bcf
bcftools view -vI samtools/bcf/${samfile[$c]%.sam}.sort.rmdup.realn.bcf > samtools/vcf/${samfile[$c]%.sam}.sort.rmdup.realn.vcf
vcfutils.pl varFilter -D1000 samtools/vcf/${samfile[$c]%.sam}.sort.rmdup.realn.vcf > samtools/vcf/${samfile[$c]%.sam}.sort.rmdup.realn.filt.vcf

java -Xmx4g -jar /opt/bio/GATK/dist/GenomeAnalysisTK.jar -R ./template/fai/${faifile%.fai} -knownSites samtools/vcf/${samfile[$c]%.sam}.sort.rmdup.realn.filt.vcf -I samtools/sort.bam/${samfile[$c]%.sam}.sort.rmdup.realn.bam -dRG ${fq_name[$c]%_?_*} -fP illumina -T CountCovariates -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate -recalFile ./output/recalibration/${samfile[$c]%.sam}.sort.rmdup.realn.csv

java -Xmx4g -jar /opt/bio/GATK/dist/GenomeAnalysisTK.jar -R ./template/fai/${faifile%.fai} -I samtools/sort.bam/${samfile[$c]%.sam}.sort.rmdup.realn.bam -dRG ${fq_name[$c]%_?_*} -fP illumina -T TableRecalibration -recalFile ./output/recalibration/${samfile[$c]%.sam}.sort.rmdup.realn.csv -o samtools/sort.bam/${samfile[$c]%.sam}.sort.rmdup.recal.bam
samtools index samtools/sort.bam/${samfile[$c]%.sam}.sort.rmdup.recal.bam
mark

echo -e "\n3. Statistics for strain ${fq_name[$c]%_?_*}\n"


reads_sum=(`wc -l ./output/sam/${paired_samfile[$c]} ./output/sam/${single_samfile[$c]}`)
reads_sum[0]=`expr ${reads_sum[4]} - 4`

reads_sum_mapped_paired=(`cut -f5 ./output/sam/${paired_samfile[$c]} | sed -e '1,2d;/^0/d' | wc -l`)
reads_sum_mapped_single=(`cut -f5 ./output/sam/${single_samfile[$c]} | sed -e '1,2d;/^0/d' | wc -l`)
reads_sum_mapped=`expr $reads_sum_mapped_paired + $reads_sum_mapped_single`

percent_mapping=`echo "scale=2;$reads_sum_mapped / $reads_sum * 100" | bc`%


#export summary files
sumfile[$c]=${samfile[$c]%sam}summary
touch ./output/statistics/${sumfile[$c]}

echo "RESULT SUMMARY FRO STRAIN ${fq_name[$c]%_?_*}: "
echo "RESULT SUMMARY FRO STRAIN ${fq_name[$c]%_?_*}:" > ./output/statistics/${sumfile[$c]}

echo "Total Reads Number= ${reads_sum[0]}"
echo "Total Reads Number= ${reads_sum[0]}" >>  ./output/statistics/${sumfile[$c]}

echo "Initially Mapped Reads Number= $reads_sum_mapped"
echo "Initially Mapped Reads Number= $reads_sum_mapped" >> ./output/statistics/${sumfile[$c]}

echo "Initially Mapped Reads Percentage= $percent_mapping"
echo "Initially Mapped Reads Percentage= $percent_mapping" >> ./output/statistics/${sumfile[$c]}


rm output/sam/${saifile1[$c]}
rm output/sam/${saifile2[$c]}
rm output/sam/${saifile_single[$c]}
rm output/sam/${single_samfile[$c]}
rm output/sam/${paired_samfile[$c]}

rm samtools/bam/${paired_samfile[$c]%.sam}.bam
rm samtools/bam/${single_samfile[$c]%.sam}.bam
rm samtools/bam/${samfile[$c]%.sam}.bam

rm samtools/sort.bam/${samfile[$c]%.sam}.sort.bam
rm samtools/sort.bam/${samfile[$c]%.sam}.sort.bam.bai
rm samtools/sort.bam/${samfile[$c]%.sam}.sort.rmdup.temp.bam
rm samtools/sort.bam/${samfile[$c]%.sam}.sort.rmdup.bam
rm samtools/sort.bam/${samfile[$c]%.sam}.sort.rmdup.bam.bai
rm samtools/sort.bam/${samfile[$c]%.sam}.sort.rmdup.metrix
rm output/recalibration/${samfile[$c]%.sam}.sort.rmdup.indel.intervals

c=`expr $c + 1`
done
