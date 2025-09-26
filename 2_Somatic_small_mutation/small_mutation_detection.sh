#!/bin/bash

family=$1


## detect mutations by Mutect2
for sample in {B,P1,P2,P3,P4,P5};do
    gatk Mutect2 -R $HOME/CHM13/chm13v2.0.fa -I $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/"$family"-"$sample".bam -I $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/control.bam -normal control --f1r2-tar-gz $HOME/mutations/"$family"/"$sample"/mutect/"$family"-"$sample"_f1r2.tar.gz -O $HOME/mutations/"$family"/"$sample"/mutect/"$family"-"$sample".vcf.gz
    gatk FilterMutectCalls -R $HOME/CHM13/chm13v2.0.fa -V $HOME/mutations/"$family"/"$sample"/mutect/"$family"-"$sample".vcf.gz -O $HOME/mutations/"$family"/"$sample"/mutect/"$family"-"$sample"_filtered.vcf.gz
done


## detect mutations by Strelka2
for sample in {B,P1,P2,P3,P4,P5};do
    mkdir -p $HOME/mutations/"$family"/"$sample"/strelka2/matched_manta
    configManta.py --normalBam $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/control.bam --tumorBam $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/"$family"-"$sample".bam  --referenceFasta $HOME/CHM13/chm13v2.0.fa --runDir $HOME/mutations/"$family"/"$sample"/strelka2/matched_manta
    $HOME/mutations/"$family"/"$sample"/strelka2/matched_manta/runWorkflow.py -m local -j 8
    configureStrelkaSomaticWorkflow.py --normalBam $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/control.bam --tumorBam $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/"$family"-"$sample".bam --referenceFasta $HOME/CHM13/chm13v2.0.fa --indelCandidates $HOME/mutations/"$family"/"$sample"/strelka2/matched_manta/results/variants/candidateSmallIndels.vcf.gz --runDir $HOME/mutations/"$family"/"$sample"/strelka2/somatic
    $HOME/mutations/"$family"/"$sample"/strelka2/somatic/runWorkflow.py -m local -j 8
done


## detect germline mutations by Strelka2 in parental genomes
for sample in {F,M};do
    configureStrelkaGermlineWorkflow.py --bam $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/"$family"-"$sample".bam  --referenceFasta $HOME/CHM13/chm13v2.0.fa --runDir $HOME/mutations/"$family"/"$sample"/germline
    $HOME/mutations/"$family"/"$sample"/germline/runWorkflow.py -m local -j 8
done

for sample in {F,M};do
    zgrep '^chr' $HOME/mutations/"$family"/"$sample"/germline/results/variants/variants.vcf.gz|awk '$7=="PASS"'|sort -V > $HOME/mutations/"$family"/"$sample"/germline/"$family"-"$sample"_germline_filter.vcf
done
