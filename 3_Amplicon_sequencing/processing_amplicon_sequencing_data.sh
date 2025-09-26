#!/bin/bash

## processing sequencing data
family=$1
tissue=$2
sample=$1"-"$2
file=$3

mkdir -p $HOME/target_PCR_data/processing_data/fastqc/"$sample"
fastqc $HOME/target_PCR_data/data/*/Data/"$file"_good_1.fq.gz $HOME/target_PCR_data/data/*/Data/"$file"_good_2.fq.gz -o $HOME/target_PCR_data/processing_data/fastqc/"$sample" -t 8

trimmomatic PE -threads 16 -phred33 $HOME/target_PCR_data/data/*/Data/"$file"_good_1.fq.gz $HOME/target_PCR_data/data/*/Data/"$file"_good_2.fq.gz $HOME/target_PCR_data/processing_data/trim/"$sample".R1.fq.gz $HOME/target_PCR_data/processing_data/trim/"$sample".R1.unpaired.fq.gz $HOME/target_PCR_data/processing_data/trim/"$sample".R2.fq.gz $HOME/target_PCR_data/processing_data/trim/"$sample".R2.unpaired.fq.gz ILLUMINACLIP:~/software/miniconda3/envs/placental_mosaicism/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:5:20 MINLEN:100

mkdir -p $HOME/target_PCR_data/processing_data/trim_fastqc/"$sample"
fastqc $HOME/target_PCR_data/processing_data/trim/"$sample".R1.fq.gz $HOME/target_PCR_data/processing_data/trim/"$sample".R2.fq.gz -o $HOME/target_PCR_data/processing_data/trim_fastqc/"$sample" -t 8

~/software/miniconda3/envs/placental_mosaicism/bin/bwa-mem2 mem -t 8 -M $HOME/CHM13/chm13v2.0.fa $HOME/target_PCR_data/processing_data/trim/"$sample".R1.fq.gz $HOME/target_PCR_data/processing_data/trim/"$sample".R2.fq.gz | samtools view -bS | samtools sort -o $HOME/target_PCR_data/processing_data/bwa/"$sample".bam

samtools view -@ 8 -q 30 -F 3340 $HOME/target_PCR_data/processing_data/bwa/"$sample".bam -b -o $HOME/target_PCR_data/processing_data/bwa_filter/"$sample".bam

samtools index $HOME/target_PCR_data/processing_data/bwa_filter/"$sample".bam > $HOME/target_PCR_data/processing_data/bwa_filter/"$sample".bam.bai


## get VAFs
awk -F '[\t_]' '{OFS="\t";print $1,$2,$2,$3,$4}' $HOME/SNV_CHM13/unmatched/"$family"/all.shared.filter.txt|awk 'length($4)==1&&length($5)==1' > $HOME/target_PCR_data/mut/"$family"/"$family".SNV.shared.bed

if [ -f $HOME/target_PCR_data/processing_data/bwa_filter/"$family""$tissue".bam ];then
    bam-readcount -w 1 $HOME/target_PCR_data/processing_data/bwa_filter/"$family""$tissue".bam -b 20 -f $HOME/CHM13/chm13v2.0.fa -l $HOME/target_PCR_data/mut/"$family"/"$family".SNV.shared.bed|awk 'BEGIN{FS=OFS="\t"}{split($6,A,":");split($7,C,":");split($8,G,":");split($9,T,":");print $1,$2,$4,A[2],C[2],G[2],T[2]}' > $HOME/target_PCR_data/mut/"$family"/all_SNV.count."$tissue".txt
    awk 'BEGIN {FS="\t"} NR==FNR {a[$1"\t"$2]=$4"\t"$5; next} ($1"\t"$2) in a {print $1"\t"$2"\t"a[$1"\t"$2]"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' $HOME/target_PCR_data/mut/"$family"/"$family".SNV.shared.bed $HOME/target_PCR_data/mut/"$family"/all_SNV.count."$tissue".txt|awk 'BEGIN {OFS="\t"} {if ($5 == 0) {print $1"_"$2"_"$3"_"$4"\t"$5} else if ($4 == "A") {print $1"_"$2"_"$3"_"$4"\t"$6/$5} else if ($4 == "C") {print $1"_"$2"_"$3"_"$4"\t"$7/$5} else if ($4 == "G") {print $1"_"$2"_"$3"_"$4"\t"$8/$5}  else if ($4 == "C") {print $1"_"$2"_"$3"_"$4"\t"$7/$5} else if ($4 == "T") {print $1"_"$2"_"$3"_"$4"\t"$9/$5}}' > $HOME/target_PCR_data/mut/"$family"/all_SNV.VAF."$tissue".txt
fi

awk -F '[\t_]' '{OFS="\t";print $1,$2,$2,$3,$4}' $HOME/SNV_CHM13/unmatched/"$family"/all.shared.filter.txt|awk '!(length($4)==1&&length($5)==1)' > $HOME/target_PCR_data/mut/"$family"/"$family".indel.shared.bed
if [ -f $HOME/target_PCR_data/processing_data/bwa_filter/"$family""$tissue".bam ];then
    python $HOME/indel_VAF.py $HOME/target_PCR_data/mut/"$family"/"$family".indel.shared.bed $HOME/target_PCR_data/processing_data/bwa_filter/"$family""$tissue".bam $HOME/target_PCR_data/mut/"$family"/all_indel.VAF."$tissue".txt $HOME/target_PCR_data/mut/"$family"/all_indel.depth."$tissue".txt
fi
