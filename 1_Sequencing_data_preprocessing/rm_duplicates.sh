#!/bin/bash

sample=$3
mkdir -p $HOME/raw_data_cutadapter_bwa_rmdup/"$sample"

samtools sort -n -@ 8 -m 2G $1 -o $HOME/raw_data_cutadapter_bwa_rmdup/"$sample"_1.bam
samtools fixmate -m -@ 8 $HOME/raw_data_cutadapter_bwa_rmdup/"$sample"_1.bam $HOME/raw_data_cutadapter_bwa_rmdup/"$sample"_2.bam
samtools sort -@ 8 -m 2G $HOME/raw_data_cutadapter_bwa_rmdup/"$sample"_2.bam -o $HOME/raw_data_cutadapter_bwa_rmdup/"$sample"_3.bam
samtools markdup -r -@ 8 $HOME/raw_data_cutadapter_bwa_rmdup/"$sample"_3.bam $2

rm $HOME/raw_data_cutadapter_bwa_rmdup/"$sample"_{1..3}.bam
rm -rf $HOME/raw_data_cutadapter_bwa_rmdup/"$sample"
