#!/bin/bash
#!/usr/bin/perl

#get localtime
ltime=`date +%Y-%m-%d-%H%M%S`

echo -n "Enter the raw FASTAQ file Directory: "
read fqd
fqd=${fqd%/}/

scoring=""

while [ "$scoring" != "yes" ] && [ "$scoring" != "no" ] && [ "$scoring" != "y" ] && [ "$scoring" != "n" ];
do
echo -n "Are they using Sanger scoring? (yes/no) "
read scoring
done

if [ "$scoring" = "y" ] || [ "$scoring" = "yes" ];
then
mkdir fastaq/fastaq_trimmed/$ltime
fi
if [ "$scoring" = "n" ] || [ "$scoring" = "no" ];
then
mkdir fastaq/sanger_scoring/$ltime
mkdir fastaq/fastaq_trimmed/$ltime
fi


#check paired-ended fq files

m=0
n=0
l=0
for fq in `ls $fqd*.fq`
do
  fq_name_ori[$m]=${fq#$fqd}
  fq_a[$m]=${fq_name_ori[$m]#*_}

  if [ "${fq_a[$m]%.*}" = "1" ];
  then
    fq1[$n]=$fq
    n=`expr $n + 1`
  fi

  if [ "${fq_a[$m]%.*}" = "2" ];
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
      if [ "${fq1[$j]%_?.*}" = "${fq2[$h]%_?.*}" ];
      then
        fq1_name[$k]=${fq1[$j]#$fqd}
        fq2_name[$k]=${fq2[$h]#$fqd}
        fq1_path[$k]=${fq1[$j]}
        fq2_path[$k]=${fq2[$h]}
        fq_name[$k]=${fq1_name[$j]%_?.*}
        k=`expr $k + 1`
        break
      else
        h=`expr $h + 1`
      fi
    done
j=`expr $j + 1`
done

echo -e "\n $k paired *.fq file(s) found! \n"



file_num=0
while [ $file_num -lt $k ]
do
echo -e "dealing with ${fq_name[$file_num]} \n"

    if [ "$scoring" = "n" ] || [ "$scoring" = "no" ];
    then
        perl ./script/sol2sanger.pl ${fq1_path[$file_num]} > fastaq/sanger_scoring/$ltime/${ltime}_${fq1_name[$file_num]%.*}_san.fq
        perl ./script/sol2sanger.pl ${fq2_path[$file_num]} > fastaq/sanger_scoring/$ltime/${ltime}_${fq2_name[$file_num]%.*}_san.fq
scythe -a /opt/bio/vsbuffalo-scythe-872a54c/truseq_adapters.fasta --quiet -q sanger -o fastaq/fastaq_trimmed/$ltime/${ltime}_${fq1_name[$file_num]%.*}_san_adaptor.fq fastaq/sanger_scoring/$ltime/${ltime}_${fq1_name[$file_num]%.*}_san.fq
scythe -a /opt/bio/vsbuffalo-scythe-872a54c/truseq_adapters.fasta --quiet -q sanger -o fastaq/fastaq_trimmed/$ltime/${ltime}_${fq2_name[$file_num]%.*}_san_adaptor.fq fastaq/sanger_scoring/$ltime/${ltime}_${fq2_name[$file_num]%.*}_san.fq

sickle pe -f fastaq/fastaq_trimmed/$ltime/${ltime}_${fq1_name[$file_num]%.*}_san_adaptor.fq -r fastaq/fastaq_trimmed/$ltime/${ltime}_${fq2_name[$file_num]%.*}_san_adaptor.fq -t sanger -o fastaq/fastaq_trimmed/$ltime/${ltime}_${fq1_name[$file_num]%.*}_san_trimmed.fq -p fastaq/fastaq_trimmed/$ltime/${ltime}_${fq2_name[$file_num]%.*}_san_trimmed.fq -s fastaq/fastaq_trimmed/$ltime/${ltime}_${fq_name[$file_num]}_single_san_trimmed.fq

rm fastaq/sanger_scoring/$ltime/${ltime}_${fq1_name[$file_num]%.*}_san.fq
rm fastaq/sanger_scoring/$ltime/${ltime}_${fq2_name[$file_num]%.*}_san.fq
rm fastaq/fastaq_trimmed/$ltime/${ltime}_${fq1_name[$file_num]%.*}_san_adaptor.fq
rm fastaq/fastaq_trimmed/$ltime/${ltime}_${fq2_name[$file_num]%.*}_san_adaptor.fq

#perl ./script/trimBWAstyle_GangMod.pl fastaq/sanger_scoring/$ltime/${ltime}_${fq%.*}_san.fq > fastaq/fastaq_trimmed/$ltime/${ltime}_${fq%.*}_san_trimmed.fq
    fi

    if [ "$scoring" = "y" ] || [ "$scoring" = "yes" ];
    then
scythe -a /opt/bio/vsbuffalo-scythe-872a54c/truseq_adapters.fasta --quiet -q sanger -o fastaq/fastaq_trimmed/$ltime/${ltime}_${fq1_name[$file_num]%.*}_san_adaptor.fq ${fq1_path[$file_num]}
scythe -a /opt/bio/vsbuffalo-scythe-872a54c/truseq_adapters.fasta --quiet -q sanger -o fastaq/fastaq_trimmed/$ltime/${ltime}_${fq2_name[$file_num]%.*}_san_adaptor.fq ${fq2_path[$file_num]}

sickle pe -f fastaq/fastaq_trimmed/$ltime/${ltime}_${fq1_name[$file_num]%.*}_san_adaptor.fq -r fastaq/fastaq_trimmed/$ltime/${ltime}_${fq2_name[$file_num]%.*}_san_adaptor.fq -t sanger -o fastaq/fastaq_trimmed/$ltime/${ltime}_${fq1_name[$file_num]%.*}_san_trimmed.fq -p fastaq/fastaq_trimmed/$ltime/${ltime}_${fq2_name[$file_num]%.*}_san_trimmed.fq -s fastaq/fastaq_trimmed/$ltime/${ltime}_${fq_name[$file_num]}_single_san_trimmed.fq

rm fastaq/fastaq_trimmed/$ltime/${ltime}_${fq1_name[$file_num]%.*}_san_adaptor.fq
rm fastaq/fastaq_trimmed/$ltime/${ltime}_${fq2_name[$file_num]%.*}_san_adaptor.fq

#perl ./script/trimBWAstyle_GangMod.pl $fq_path > fastaq/fastaq_trimmed/$ltime/${ltime}_${fq%.*}_san_trimmed.fq
    fi

file_num=`expr $file_num + 1`
done
