#!/bin/bash

family=$1
mkdir -p $HOME/mutations/"$family"/filter


## find clustered mutations in unfiltered mutation sets
for sample in {B,P1,P2,P3,P4,P5};do
        cat <(cat $HOME/mutations/"$family"/"$sample"/somaticseq/"$family"-"$sample".strelka2.snv.vcf|grep -v "#"|cut -f 1,2,4,5) <(cat $HOME/mutations/"$family"/"$sample"/somaticseq/"$family"-"$sample".strelka2.indel.vcf|grep -v "#"|cut -f 1,2,4,5) <(cat $HOME/mutations/"$family"/"$sample"/somaticseq/"$family"-"$sample".mutect2.vcf|grep -v "#"|cut -f 1,2,4,5)|sort -V|uniq|awk -F '[\t_]' '{chr = $1;pos = $2;if (NR == 1) {prev_chr = chr;prev_pos = pos;prev_line = $0;next} if (chr == prev_chr && (pos - prev_pos) <= 150) {print prev_line;print $0} prev_chr = chr;prev_pos = pos;prev_line = $0}'|uniq > $HOME/mutations/"$family"/filter/"$family"."$sample".filter.txt
done

cat $HOME/mutations/"$family"/filter/"$family".{B,P1,P2,P3,P4,P5}.filter.txt|sort -V|uniq|awk '{print $1"_"$2"\t"$3"\t"$4}' > $HOME/mutations/"$family"/filter/"$family".filter.site.txt


## find clustered mutations in filtered mutation sets
cat $HOME/mutations/"$family"/all.shared.txt $HOME/mutations/"$family"/all.unique.txt | sort -V | awk -F '[\t_]' '{print $1"_"$2}' |while read line;do grep -w $line -A2 -B2 $HOME/mutations/"$family"/filter/"$family".filter.site.txt | awk -F '[\t_]' '{chr = $1;pos = $2;if (NR == 1) {prev_chr = chr;prev_pos = pos;prev_line = $0;next} if (chr == prev_chr && (pos - prev_pos) <= 150) {print prev_line;print $0} prev_chr = chr;prev_pos = pos;prev_line = $0}';done|sort -V|uniq|awk -F '[\t_]' '{print $1"\t"$2"\t"$2"\t"$3"\t"$4}' > $HOME/mutations/"$family"/filter/"$family".filter.test.txt


## find germline SNVs
awk 'length($4)==1 && length($5)==1' $HOME/mutations/"$family"/filter/"$family".filter.test.txt > $HOME/mutations/"$family"/filter/"$family".filter.test.SNV.txt

for sample in {M,F,B,P1,P2,P3,P4,P5};do
    bam-readcount -w 1 $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/"$family"-"$sample".bam -b 20 -f $HOME/CHM13/chm13v2.0.fa -l $HOME/mutations/"$family"/filter/"$family".filter.test.SNV.txt|awk 'BEGIN{FS=OFS="\t"}{split($6,A,":");split($7,C,":");split($8,G,":");split($9,T,":");print $1,$2,$4,A[2],C[2],G[2],T[2]}' > $HOME/mutations/"$family"/filter/"$family"-"$sample".SNV.count.txt
done

for sample in {F,M};do
    awk 'BEGIN {FS="\t"} NR==FNR {a[$1"\t"$2]=$4"\t"$5; next} ($1"\t"$2) in a {print $1"\t"$2"\t"a[$1"\t"$2]"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' $HOME/mutations/"$family"/filter/"$family".filter.test.SNV.txt $HOME/mutations/"$family"/filter/"$family"-"$sample".SNV.count.txt|awk 'BEGIN {OFS="\t"} {if (($4 == "A" && $6!=0) || ($4 == "C" && $7!=0) || ($4 == "G" && $8!=0) || ($4 == "T" && $9!=0)) {$NF=$NF"\t""T"} else {$NF=$NF"\t""F"} print }' > $HOME/mutations/"$family"/filter/"$family"-"$sample".SNV.validation.txt
done

for sample in {B,P1,P2,P3,P4,P5};do
    awk 'BEGIN {FS="\t"} NR==FNR {a[$1"\t"$2]=$4"\t"$5; next} ($1"\t"$2) in a {print $1"\t"$2"\t"a[$1"\t"$2]"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' $HOME/mutations/"$family"/filter/"$family".filter.test.SNV.txt $HOME/mutations/"$family"/filter/"$family"-"$sample".SNV.count.txt|awk 'BEGIN {OFS="\t"} {if ($5!=0 && (($4 == "A" && $6/$5 > 0.25) || ($4 == "C" && $7/$5 > 0.25) || ($4 == "G" && $8/$5 > 0.25) || ($4 == "T" && $9/$5 > 0.25))) {$NF=$NF"\t""T"} else {$NF=$NF"\t""F"} print }' > $HOME/mutations/"$family"/filter/"$family"-"$sample".SNV.validation.txt
done

cd $HOME/mutations/"$family"/filter
awk '{if(FNR==NR){s[$1"_"$2]++;val[$1"_"$2]=$10}else{if(s[$1"_"$2]==1){print $0,val[$1"_"$2]}}}' "$family"-P5.SNV.validation.txt <(awk '{if(FNR==NR){s[$1"_"$2]++;val[$1"_"$2]=$10}else{if(s[$1"_"$2]==1){print $0,val[$1"_"$2]}}}' "$family"-P4.SNV.validation.txt <(awk '{if(FNR==NR){s[$1"_"$2]++;val[$1"_"$2]=$10}else{if(s[$1"_"$2]==1){print $0,val[$1"_"$2]}}}' "$family"-P3.SNV.validation.txt <(awk '{if(FNR==NR){s[$1"_"$2]++;val[$1"_"$2]=$10}else{if(s[$1"_"$2]==1){print $0,val[$1"_"$2]}}}' "$family"-P2.SNV.validation.txt <(awk '{if(FNR==NR){s[$1"_"$2]++;val[$1"_"$2]=$10}else{if(s[$1"_"$2]==1){print $0,val[$1"_"$2]}}}' "$family"-P1.SNV.validation.txt <(awk '{if(FNR==NR){s[$1"_"$2]++;val[$1"_"$2]=$10}else{if(s[$1"_"$2]==1){print $0,val[$1"_"$2]}}}' "$family"-B.SNV.validation.txt <(awk '{if(FNR==NR){s[$1"_"$2]++;val[$1"_"$2]=$10}else{if(s[$1"_"$2]==1){print $1,$2,$3,$4,$10,val[$1"_"$2]}}}' "$family"-F.SNV.validation.txt "$family"-M.SNV.validation.txt)))))) > "$family"-all.SNV.validation.txt


## find germline indels
awk 'length($4)!=1 || length($5)!=1' $HOME/mutations/"$family"/filter/"$family".filter.test.txt > $HOME/mutations/"$family"/filter/"$family".filter.test.indel.txt

python indel_VAF.py $HOME/mutations/"$family"/filter/"$family".filter.test.indel.txt $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/"$family"-M.bam $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/"$family"-F.bam $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/"$family"-B.bam $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/"$family"-P1.bam $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/"$family"-P2.bam $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/"$family"-P3.bam $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/"$family"-P4.bam $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/"$family"-P5.bam $HOME/mutations/"$family"/filter/"$family"-all.indel.validation.txt


## find candidate clustered mutations
cat <(awk '($5=="T" || $6=="T") && ($7=="T" && $8=="T" && $9=="T" && $10=="T" && $11=="T" && $12=="T")' $HOME/mutations/"$family"/filter/"$family"-all.SNV.validation.txt|awk 'length($3)==1 && length($4)==1'|cut -f 1,2,3,4) <(awk '($2>0 || $3>0) && ($4>0.25 && $5>0.25 && $6>0.25 && $7>0.25 && $8>0.25 && $9>0.25)' $HOME/mutations/"$family"/filter/"$family"-all.indel.validation.txt|awk -F '[\t_]' '{OFS="\t";print $1,$2,$3,$4}')|sort -V > $HOME/mutations/"$family"/filter/"$family"-all.germline.txt

awk '{if(FNR==NR){s[$1"_"$2"_"$3"_"$4]++}else{if(s[$1"_"$2"_"$4"_"$5]==0){print $0}}}' $HOME/mutations/"$family"/filter/"$family"-all.germline.txt $HOME/mutations/"$family"/filter/"$family".filter.test.txt|awk -F '[\t_]' '{chr = $1;pos = $2;if (NR == 1) {prev_chr = chr;prev_pos = pos;prev_line = $0;next} if (chr == prev_chr && (pos - prev_pos) <= 150) {print prev_line;print $0} prev_chr = chr;prev_pos = pos;prev_line = $0}'|sort -V|uniq > $HOME/mutations/"$family"/filter/"$family".filter.txt


## generate the final mutation sets
cat <(ls $HOME/mutations/check_SNV/check_SNV_family/{M1,M2,M3,M4,M5,M6,N1,N2,N3,C1,C2,C3}.validation.txt|grep -v "$family"|xargs cat|awk -F "\t" '{count=0; for(i=2;i<=NF;i++) if($i=="F") count++; if(count==NF-1) print $1}'|sort -V|uniq -c|awk '$1==11{print $2}') <(ls $HOME/mutations/check_indel/check_indel_family/{M1,M2,M3,M4,M5,M6,N1,N2,N3,C1,C2,C3}.check_indel.txt|grep -v "$family"|xargs cat|awk -F "\t" '{count=0; for(i=5;i<=NF;i++) if($i=="F") count++; if(count==NF-4) print $1"_"$2"_"$3"_"$4}'|sort -V|uniq -c|awk '$1==11{print $2}') > $HOME/mutations/"$family"/new.mut.txt
grep -wFf $HOME/mutations/"$family"/new.mut.txt <(awk -F '[\t_]' '{if(FNR==NR){s[$1"_"$2]++}else{if(s[$1"_"$2]==0){print $0}}}' $HOME/mutations/"$family"/filter/"$family".filter.txt $HOME/mutations/"$family"/all.shared.txt)|sort -V|awk -F '[\t_]' '{chr[NR] = $1;positions[NR] = $2;lines[NR] = $0}END {for (i = 1; i <= NR; i++) {if ((i == 1 && (chr[i] != chr[i+1] || positions[i+1] - positions[i] > 150)) || (i > 1 && (chr[i] != chr[i-1] || positions[i] - positions[i-1] > 150) && (i == NR || chr[i] != chr[i+1] || positions[i+1] - positions[i] > 150))) {print lines[i]}}}'|sort -V|uniq > $HOME/mutations/"$family"/all.shared.filter.txt
grep -wFf $HOME/mutations/"$family"/new.mut.txt <(awk -F '[\t_]' '{if(FNR==NR){s[$1"_"$2]++}else{if(s[$1"_"$2]==0){print $0}}}' $HOME/mutations/"$family"/filter/"$family".filter.txt $HOME/mutations/"$family"/all.unique.txt)|sort -V|awk -F '[\t_]' '{chr[NR] = $1;positions[NR] = $2;lines[NR] = $0}END {for (i = 1; i <= NR; i++) {if ((i == 1 && (chr[i] != chr[i+1] || positions[i+1] - positions[i] > 150)) || (i > 1 && (chr[i] != chr[i-1] || positions[i] - positions[i-1] > 150) && (i == NR || chr[i] != chr[i+1] || positions[i+1] - positions[i] > 150))) {print lines[i]}}}'|sort -V|uniq > $HOME/mutations/"$family"/all.unique.filter.txt
