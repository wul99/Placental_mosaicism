#!/bin/bash

family=$1


## find SNVs present in more than one family
mkdir -p $HOME/mutations/check_SNV/check_SNV_family

cat $HOME/mutations/*/all.unique.txt $HOME/mutations/*/all.shared.txt|awk -F '[\t_]' 'length($3)==1&&length($4)==1'|sort -V|awk -F '[\t_]' '{OFS="\t";print $1,$2,$2,$3,$4}' > $HOME/mutations/check_SNV/all.SNV.bed

for sample in {M,F,B,P1,P2,P3,P4,P5};do
    bam-readcount -w 1 $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/"$family"-"$sample".bam -b 20 -f $HOME/CHM13/chm13v2.0.fa -l $HOME/mutations/check_SNV/all.SNV.bed|awk 'BEGIN{FS=OFS="\t"}{split($6,A,":");split($7,C,":");split($8,G,":");split($9,T,":");print $1,$2,A[2],C[2],G[2],T[2]}' > $HOME/mutations/check_SNV/check_SNV_family/"$family"-"$sample".txt
    awk 'BEGIN {FS="\t"} NR==FNR {a[$1"\t"$2]=$4"\t"$5; next} ($1"\t"$2) in a {print $1"\t"$2"\t"a[$1"\t"$2]"\t"$3"\t"$4"\t"$5"\t"$6}' $HOME/mutations/check_SNV/all.SNV.bed $HOME/mutations/check_SNV/check_SNV_family/"$family"-"$sample".txt|awk 'BEGIN {OFS="\t"} {if (($4 == "A" && $5!=0) || ($4 == "C" && $6!=0) || ($4 == "G" && $7!=0) || ($4 == "T" && $8!=0)) {print $1"_"$2"_"$3"_"$4"\t""T"} else {print $1"_"$2"_"$3"_"$4"\t""F"}}' > $HOME/mutations/check_SNV/check_SNV_family/"$family"-"$sample".validation.txt
done
cd $HOME/mutations/check_SNV/check_SNV_family
awk 'BEGIN {FS="\t"} NR==FNR {a[$1]=$2; next} $1 in a {print $0"\t"a[$1]}' "$family"-P5.validation.txt <(awk 'BEGIN {FS="\t"} NR==FNR {a[$1]=$2; next} $1 in a {print $0"\t"a[$1]}' "$family"-P4.validation.txt <(awk 'BEGIN {FS="\t"} NR==FNR {a[$1]=$2; next} $1 in a {print $0"\t"a[$1]}' "$family"-P3.validation.txt <(awk 'BEGIN {FS="\t"} NR==FNR {a[$1]=$2; next} $1 in a {print $0"\t"a[$1]}' "$family"-P2.validation.txt <(awk 'BEGIN {FS="\t"} NR==FNR {a[$1]=$2; next} $1 in a {print $0"\t"a[$1]}' "$family"-P1.validation.txt <(awk 'BEGIN {FS="\t"} NR==FNR {a[$1]=$2; next} $1 in a {print $0"\t"a[$1]}' "$family"-B.validation.txt <(awk 'BEGIN {FS="\t"} NR==FNR {a[$1]=$2; next} $1 in a {print $0"\t"a[$1]}' "$family"-F.validation.txt <(awk 'BEGIN {FS="\t"} NR==FNR {a[$1]=$2; next} $1 in a {print $0"\t"a[$1]}' "$family"-M.validation.txt <(cat "$family"-*.validation.txt|cut -f1|sort -V|uniq)))))))) > "$family".validation.txt


## find indels present in more than one family
mkdir -p $HOME/mutations/check_indel/check_indel_family

cat $HOME/mutations/*/all.unique.txt $HOME/mutations/*/all.shared.txt|awk -F '[\t_]' '!(length($3)==1&&length($4)==1)'|cut -f1 |awk -F '[\t_]' '{OFS="\t";print $1,$2,".",$3,$4}' > $HOME/mutations/check_indel/all.indel.txt

python check_indel.py $HOME/mutations/check_indel/all.indel.txt $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/"$family"-M.bam $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/"$family"-F.bam $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/"$family"-B.bam $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/"$family"-P1.bam $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/"$family"-P2.bam $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/"$family"-P3.bam $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/"$family"-P4.bam $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/"$family"-P5.bam $HOME/mutations/check_indel/check_indel_family/"$family".check_indel.txt
