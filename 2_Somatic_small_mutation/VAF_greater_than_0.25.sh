#!/bin/bash

family=$1
mkdir -p $HOME/mutations/"$family"/new_VAF


## calculate SNV VAFs
cat $HOME/mutations/"$family"/*/somaticseq/"$family"-*.sSNV.MF_filter.vcf|awk '{OFS="\t";print $1,$2,$2,$4,$5}'|sort -V|uniq > $HOME/mutations/"$family"/new_VAF/all_SNV.bed

for sample in {B,P1,P2,P3,P4,P5};do
    bam-readcount -w 1 $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/"$family"-"$sample".bam -b 20 -f $HOME/CHM13/chm13v2.0.fa -l $HOME/mutations/"$family"/new_VAF/all_SNV.bed|awk 'BEGIN{FS=OFS="\t"}{split($6,A,":");split($7,C,":");split($8,G,":");split($9,T,":");print $1,$2,$4,A[2],C[2],G[2],T[2]}' > $HOME/mutations/"$family"/new_VAF/all_SNV.count."$sample".txt
    awk 'BEGIN {FS="\t"} NR==FNR {a[$1"\t"$2]=$4"\t"$5; next} ($1"\t"$2) in a {print $1"\t"$2"\t"a[$1"\t"$2]"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' $HOME/mutations/"$family"/new_VAF/all_SNV.bed $HOME/mutations/"$family"/new_VAF/all_SNV.count."$sample".txt|awk 'BEGIN {OFS="\t"} {if ($5 == 0) {print $1"_"$2"_"$3"_"$4"\t"$5} else if ($4 == "A") {print $1"_"$2"_"$3"_"$4"\t"$6/$5} else if ($4 == "C") {print $1"_"$2"_"$3"_"$4"\t"$7/$5} else if ($4 == "G") {print $1"_"$2"_"$3"_"$4"\t"$8/$5}  else if ($4 == "C") {print $1"_"$2"_"$3"_"$4"\t"$7/$5} else if ($4 == "T") {print $1"_"$2"_"$3"_"$4"\t"$9/$5}}' > $HOME/mutations/"$family"/new_VAF/all_SNV.VAF."$sample".txt
done

cd $HOME/mutations/"$family"/new_VAF
cat $HOME/mutations/"$family"/new_VAF/all_SNV.VAF.{B,P1,P2,P3,P4,P5}.txt|cut -f1|sort -V|uniq|awk 'BEGIN{while(getline<"all_SNV.VAF.B.txt")a[$1]=$2}{OFS="\t";print $0,a[$1]}'|awk 'BEGIN{while(getline<"all_SNV.VAF.P1.txt")a[$1]=$2}{OFS="\t";print $0,a[$1]}'|awk 'BEGIN{while(getline<"all_SNV.VAF.P2.txt")a[$1]=$2}{OFS="\t";print $0,a[$1]}'|awk 'BEGIN{while(getline<"all_SNV.VAF.P3.txt")a[$1]=$2}{OFS="\t";print $0,a[$1]}'|awk 'BEGIN{while(getline<"all_SNV.VAF.P4.txt")a[$1]=$2}{OFS="\t";print $0,a[$1]}'|awk 'BEGIN{while(getline<"all_SNV.VAF.P5.txt")a[$1]=$2}{OFS="\t";print $0,a[$1]}' > $HOME/mutations/"$family"/new_VAF/all_SNV.VAF.txt

awk 'FNR==NR{a[$1,$2]=$3; next} {print $0"\t"(a[$1,$2] ? a[$1,$2] : 0)}' $HOME/mutations/"$family"/new_VAF/all_SNV.count.P5.txt <(awk 'FNR==NR{a[$1,$2]=$3; next} {print $0"\t"(a[$1,$2] ? a[$1,$2] : 0)}' $HOME/mutations/"$family"/new_VAF/all_SNV.count.P4.txt <(awk 'FNR==NR{a[$1,$2]=$3; next} {print $0"\t"(a[$1,$2] ? a[$1,$2] : 0)}' $HOME/mutations/"$family"/new_VAF/all_SNV.count.P3.txt <(awk 'FNR==NR{a[$1,$2]=$3; next} {print $0"\t"(a[$1,$2] ? a[$1,$2] : 0)}' $HOME/mutations/"$family"/new_VAF/all_SNV.count.P2.txt <(awk 'FNR==NR{a[$1,$2]=$3; next} {print $0"\t"(a[$1,$2] ? a[$1,$2] : 0)}' $HOME/mutations/"$family"/new_VAF/all_SNV.count.P1.txt <(awk 'FNR==NR{a[$1,$2]=$3; next} {print $0"\t"(a[$1,$2] ? a[$1,$2] : 0)}' $HOME/mutations/"$family"/new_VAF/all_SNV.count.B.txt $HOME/mutations/"$family"/new_VAF/all_SNV.bed)))))|awk '{OFS="\t";print $1"_"$2"_"$4"_"$5,$6,$7,$8,$9,$10,$11}' > $HOME/mutations/"$family"/new_VAF/all_SNV.depth.txt


## calculate indel VAFs
cat $HOME/mutations/"$family"/*/somaticseq/"$family"-*.sINDEL.MF_filter.vcf|awk '{OFS="\t";print $1,$2,$2,$4,$5}'|sort -V|uniq > $HOME/mutations/"$family"/new_VAF/all_indel.bed

python indel_VAF.py $HOME/mutations/"$family"/new_VAF/all_indel.bed $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/"$family"-B.bam $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/"$family"-P1.bam $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/"$family"-P2.bam $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/"$family"-P3.bam $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/"$family"-P4.bam $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/"$family"-P5.bam $HOME/mutations/"$family"/new_VAF/all_indel.VAF.txt $HOME/mutations/"$family"/new_VAF/all_indel.depth.txt


## filter somatic mutations by VAFs
cat <(awk '{if(FNR==NR){s[$1]++}else{if(s[$1]==1){print $0}}}' <(awk '{count=0; for(i=2;i<=NF;i++) if($i >= 20 && $i <= 120) count++; if(count == NF-1) print}' $HOME/mutations/"$family"/new_VAF/all_indel.depth.txt) <(awk '{count=0; for(i=2;i<=NF;i++) if($i > 0) count++; if(count >= 2) print}' $HOME/mutations/"$family"/new_VAF/all_indel.VAF.txt)|awk '{count=0; for(i=2;i<=NF;i++) if($i >= 0.05) count++; if(count > 0) print}') <(awk '{if(FNR==NR){s[$1]++}else{if(s[$1]==1){print $0}}}' <(awk '{count=0; for(i=2;i<=NF;i++) if($i >= 20 && $i <= 120) count++; if(count == NF-1) print}' $HOME/mutations/"$family"/new_VAF/all_SNV.depth.txt) <(awk '{count=0; for(i=2;i<=NF;i++) if($i > 0) count++; if(count >= 2) print}' $HOME/mutations/"$family"/new_VAF/all_SNV.VAF.txt)|awk '{count=0; for(i=2;i<=NF;i++) if($i >= 0.05) count++; if(count > 0) print}')|awk '{count=0; for(i=2;i<=NF;i++) if($i > 0.25) count++; if(count != NF-1) print}'|sort -V > $HOME/mutations/"$family"/all.shared.txt
cat <(awk '{if(FNR==NR){s[$1]++}else{if(s[$1]==1){print $0}}}' <(awk '{count=0; for(i=2;i<=NF;i++) if($i >= 20 && $i <= 120) count++; if(count == NF-1) print}' $HOME/mutations/"$family"/new_VAF/all_indel.depth.txt) <(awk '{count=0; for(i=2;i<=NF;i++) if($i > 0) count++; if(count == 1) print}' $HOME/mutations/"$family"/new_VAF/all_indel.VAF.txt)) <(awk '{if(FNR==NR){s[$1]++}else{if(s[$1]==1){print $0}}}' <(awk '{count=0; for(i=2;i<=NF;i++) if($i >= 20 && $i <= 120) count++; if(count == NF-1) print}' $HOME/mutations/"$family"/new_VAF/all_SNV.depth.txt) <(awk '{count=0; for(i=2;i<=NF;i++) if($i > 0) count++; if(count == 1) print}' $HOME/mutations/"$family"/new_VAF/all_SNV.VAF.txt))|awk '{count=0; for(i=2;i<=NF;i++) if($i >= 0.05) count++; if(count > 0) print}'|sort -V > $HOME/mutations/"$family"/all.unique.txt
