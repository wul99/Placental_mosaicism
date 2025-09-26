#!/bin/bash

family=$1


## detect mutations by GRIDSS
for sample in {B,P1,P2,P3,P4,P5};do
    gridss -r $HOME/SV/gridss/ref/chm13v2.0.fa -w $HOME/SV/gridss/tmp -j ~/software/gridss-master/gridss-2.13.2-gridss-jar-with-dependencies.jar -t 8 -o $HOME/SV/gridss/"$family"/"$sample"/matched_gridss/"$family"-"$sample".vcf -b $HOME/SV/gridss/exclude_list.bed $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/control.bam $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/"$family"-"$sample".bam
    Rscript simple-event-annotation.R $HOME/SV/gridss/"$family"/"$sample"/"$family"-"$sample".vcf.gz $HOME/SV/gridss/"$family"/"$sample"/"$family"-"$sample"_annotated.vcf $HOME/SV/gridss/"$family"/"$sample"/"$family"-"$sample"_filter.bed
done

gridss --reference $HOME/SV/gridss/ref/chm13v2.0.fa --output $HOME/SV/gridss/$family/M/$family-M.vcf.gz --threads 8 --jar ~/software/gridss-master/gridss-2.13.2-gridss-jar-with-dependencies.jar --workingdir $HOME/SV/gridss/tmp $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/$family-M.bam
Rscript $HOME/SV/gridss/simple-event-annotation.MF.R $HOME/SV/gridss/"$family"/M/"$family"-M.vcf.gz $HOME/SV/gridss/"$family"/M/"$family"-M_annotated.vcf $HOME/SV/gridss/"$family"/M/"$family"-M_filter.bed

gridss --reference $HOME/SV/gridss/ref/chm13v2.0.fa --output $HOME/SV/gridss/"$family"/F/"$family"-F.vcf.gz --threads 8 --jar ~/software/gridss-master/gridss-2.13.2-gridss-jar-with-dependencies.jar --workingdir $HOME/SV/gridss/tmp $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/"$family"-F.bam
Rscript $HOME/SV/gridss/simple-event-annotation.MF.R $HOME/SV/gridss/"$family"/F/"$family"-F.vcf.gz $HOME/SV/gridss/"$family"/F/"$family"-F_annotated.vcf $HOME/SV/gridss/"$family"/F/"$family"-F_filter.bed


## detect mutations by Manta
for sample in {B,P1,P2,P3,P4,P5};do
    configManta.py --normalBam $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/control.bam --tumorBam $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/"$family"-"$sample".bam  --referenceFasta $HOME/CHM13/chm13v2.0.fa --runDir $HOME/SV/manta/"$family"/"$sample"/matched_manta
    $HOME/SV/manta/"$family"/"$sample"/matched_manta/runWorkflow.py -m local -j 8
done

configManta.py --bam $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/"$family"-F.bam --referenceFasta $HOME/CHM13/chm13v2.0.fa --runDir $HOME/SV/manta/"$family"/F/manta
$HOME/SV/manta/"$family"/F/manta/runWorkflow.py -m local -j 8

configManta.py --bam $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/"$family"-M.bam --referenceFasta $HOME/CHM13/chm13v2.0.fa --runDir $HOME/SV/manta/"$family"/M/manta
$HOME/SV/manta/"$family"/M/manta/runWorkflow.py -m local -j 8
