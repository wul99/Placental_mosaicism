#!/bin/bash

family=$1
sample=$2


## filter germline mutations detected in parental genomes by Strelka2
zcat $HOME/mutations/"$family"/"$sample"/strelka2/somatic/results/variants/somatic.*.vcf.gz|grep '^chr'|awk '$7=="PASS"'|sort -V > $HOME/mutations/"$family"/"$sample"/strelka2/"$family"-"$sample"_somatic_filter.vcf
awk -F'\t' '{if(FNR==NR){s[$1"_"$2"_"$4"_"$5]++}else{if(s[$1"_"$2"_"$4"_"$5]==0){print $0}}}' $HOME/mutations/"$family"/"$sample"/germline/"$family"-M_germline_filter.vcf $HOME/mutations/"$family"/"$sample"/strelka2/"$family"-"$sample"_somatic_filter.vcf > $HOME/mutations/"$family"/"$sample"/strelka2/"$family"-"$sample"_somatic_filter_M.vcf
awk -F'\t' '{if(FNR==NR){s[$1"_"$2"_"$4"_"$5]++}else{if(s[$1"_"$2"_"$4"_"$5]==0){print $0}}}' $HOME/mutations/"$family"/"$sample"/germline/"$family"-F_germline_filter.vcf $HOME/mutations/"$family"/"$sample"/strelka2/"$family"-"$sample"_somatic_filter_M.vcf > $HOME/mutations/"$family"/"$sample"/strelka2/"$family"-"$sample"_somatic_filter_MF.vcf


## integrate the mutation sets from both algorithms and apply additional quality filtering by SomaticSeq
cat <(zgrep "#" $HOME/mutations/"$family"/"$sample"/strelka2/somatic/results/variants/somatic.snvs.vcf.gz) <(grep "DP=" $HOME/mutations/"$family"/"$sample"/strelka2/"$family"-"$sample"_somatic_filter_MF.vcf) > $HOME/mutations/"$family"/"$sample"/somaticseq/"$family"-"$sample".strelka2.snv.vcf
cat <(zgrep "#" $HOME/mutations/"$family"/"$sample"/strelka2/somatic/results/variants/somatic.indels.vcf.gz) <(grep -v "DP=" $HOME/mutations/"$family"/"$sample"/strelka2/"$family"-"$sample"_somatic_filter_MF.vcf) > $HOME/mutations/"$family"/"$sample"/somaticseq/"$family"-"$sample".strelka2.indel.vcf
cat <(zgrep "#" $HOME/mutations/"$family"/"$sample"/mutect/"$family"-"$sample"_filtered.vcf.gz) <(zgrep -v "#" $HOME/mutations/"$family"/"$sample"/mutect/"$family"-"$sample"_filtered.vcf.gz|awk '$7=="PASS"') > $HOME/mutations/"$family"/"$sample"/somaticseq/"$family"-"$sample".mutect2.vcf

somaticseq_parallel.py --output-directory $HOME/mutations/"$family"/"$sample"/somaticseq --genome-reference $HOME/CHM13/chm13v2.0.fa --minimum-num-callers 2 --threads 6 --minimum-mapping-quality 30 --inclusion-region $HOME/all_Regions.bed --minimum-base-quality 20 paired --tumor-bam-file $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/"$family"-"$sample".bam --normal-bam-file $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/control.bam --mutect2-vcf $HOME/mutations/"$family"/"$sample"/somaticseq/"$family"-"$sample".mutect2.vcf --strelka-snv $HOME/mutations/"$family"/"$sample"/somaticseq/"$family"-"$sample".strelka2.snv.vcf --strelka-indel $HOME/mutations/"$family"/"$sample"/somaticseq/"$family"-"$sample".strelka2.indel.vcf


## find candidate germline mutations in detected somatic mutations
awk -F "[,:]" '$37+$38+$39+$40>=20' $HOME/mutations/"$family"/"$sample"/somaticseq/Consensus.sSNV.vcf|awk '$7=="PASS"{OFS="\t";print $1,$2,$2}' > $HOME/mutations/"$family"/"$sample"/somaticseq/"$family"-"$sample".pass.snv.bed

bam-readcount -w 1 $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/"$family"-M.bam -b 20 -f $HOME/CHM13/chm13v2.0.fa -l $HOME/mutations/"$family"/"$sample"/somaticseq/"$family"-"$sample".pass.snv.bed|awk 'BEGIN{FS=OFS="\t"}{split($6,A,":");split($7,C,":");split($8,G,":");split($9,T,":");print $1,$2,$4,A[2],C[2],G[2],T[2]}' > $HOME/mutations/"$family"/"$sample"/somaticseq/"$family"-M.pass.snv.Illumina_count.txt
awk 'BEGIN {FS="\t"} NR==FNR {a[$1"\t"$2]=$4"\t"$5; next} ($1"\t"$2) in a {print $1"\t"$2"\t"a[$1"\t"$2]"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' $HOME/mutations/"$family"/"$sample"/somaticseq/Consensus.sSNV.vcf $HOME/mutations/"$family"/"$sample"/somaticseq/"$family"-M.pass.snv.Illumina_count.txt|awk 'BEGIN {OFS="\t"} {if (($4 == "A" && $6==0) || ($4 == "C" && $7==0) || ($4 == "G" && $8==0) || ($4 == "T" && $9==0)) {$NF=$NF"\t""T"} else {$NF=$NF"\t""F"} print }' > $HOME/mutations/"$family"/"$sample"/somaticseq/"$family"-M.pass.snv.Illumina_count.type.validation.txt

bam-readcount -w 1 $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/"$family"-F.bam -b 20 -f $HOME/CHM13/chm13v2.0.fa -l $HOME/mutations/"$family"/"$sample"/somaticseq/"$family"-"$sample".pass.snv.bed|awk 'BEGIN{FS=OFS="\t"}{split($6,A,":");split($7,C,":");split($8,G,":");split($9,T,":");print $1,$2,$4,A[2],C[2],G[2],T[2]}' > $HOME/mutations/"$family"/"$sample"/somaticseq/"$family"-F.pass.snv.Illumina_count.txt
awk 'BEGIN {FS="\t"} NR==FNR {a[$1"\t"$2]=$4"\t"$5; next} ($1"\t"$2) in a {print $1"\t"$2"\t"a[$1"\t"$2]"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' $HOME/mutations/"$family"/"$sample"/somaticseq/Consensus.sSNV.vcf $HOME/mutations/"$family"/"$sample"/somaticseq/"$family"-F.pass.snv.Illumina_count.txt|awk 'BEGIN {OFS="\t"} {if (($4 == "A" && $6==0) || ($4 == "C" && $7==0) || ($4 == "G" && $8==0) || ($4 == "T" && $9==0)) {$NF=$NF"\t""T"} else {$NF=$NF"\t""F"} print }' > $HOME/mutations/"$family"/"$sample"/somaticseq/"$family"-F.pass.snv.Illumina_count.type.validation.txt

awk -F "[,:]" '$37+$38+$39+$40>=20' $HOME/mutations/"$family"/"$sample"/somaticseq/Consensus.sINDEL.vcf > $HOME/mutations/"$family"/"$sample"/somaticseq/sINDEL.vcf
python check_indel.py $HOME/mutations/"$family"/"$sample"/somaticseq/sINDEL.vcf $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/"$family"-M.bam $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/"$family"-F.bam $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/"$family"-"$sample".bam $HOME/mutations/"$family"/"$sample"/somaticseq/check_indel_MF.txt
awk '{if(FNR==NR){s[$1"_"$2"_"$3"_"$4]++}else{if(s[$1"_"$2"_"$4"_"$5]==1){print $0}}}' <(awk '$5=="F" && $6=="F" && $7=="T"' $HOME/mutations/"$family"/"$sample"/somaticseq/check_indel_MF.txt) $HOME/mutations/"$family"/"$sample"/somaticseq/Consensus.sINDEL.vcf|awk -F "[,:]" '$37+$38+$39+$40>=20'|awk -F ":" '$NF>0'|awk '{print $1":"$2"-"$2}' > $HOME/mutations/"$family"/"$sample"/somaticseq/"$family"-"$sample".sINDEL.txt
cat $HOME/mutations/"$family"/"$sample"/somaticseq/"$family"-"$sample".sINDEL.txt | while read line;do samtools depth -r $line $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/"$family"-F.bam $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/"$family"-M.bam $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/"$family"-"$sample".bam;done > $HOME/mutations/"$family"/"$sample"/somaticseq/"$family"-"$sample".sINDEL.depth.txt


## filter somatic mutations
awk '{if(FNR==NR){s[$1"_"$2"_"$3"_"$4]++}else{if(s[$1"_"$2"_"$4"_"$5]==1){print $0}}}' <(cat <(awk '$10=="T" && $5>=10' $HOME/mutations/"$family"/"$sample"/somaticseq/"$family"-F.pass.snv.Illumina_count.type.validation.txt|cut -f 1,2,3,4) <(awk '$10=="T" && $5>=10' $HOME/mutations/"$family"/"$sample"/somaticseq/"$family"-M.pass.snv.Illumina_count.type.validation.txt|cut -f 1,2,3,4)|sort|uniq -c|awk '$1==2{OFS="\t";print $2,$3,$4,$5}') $HOME/mutations/"$family"/"$sample"/somaticseq/Consensus.sSNV.vcf|awk -F "[,:]" '$37+$38+$39+$40>=20' > $HOME/mutations/"$family"/"$sample"/somaticseq/"$family"-"$sample".sSNV.MF_filter.vcf
awk '{if(FNR==NR){s[$1"_"$2]++}else{if(s[$1"_"$2]==1){print $0}}}' <(awk '$3>=10 && $4>=10 && $5>=20' $HOME/mutations/"$family"/"$sample"/somaticseq/"$family"-"$sample".sINDEL.depth.txt) $HOME/mutations/"$family"/"$sample"/somaticseq/Consensus.sINDEL.vcf > $HOME/mutations/"$family"/"$sample"/somaticseq/"$family"-"$sample".sINDEL.MF_filter.vcf
