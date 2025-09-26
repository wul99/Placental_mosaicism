#!/bin/bash

## adapters and low-quality bases were trimmed using Trimmomatic
cat $HOME/gene_expression/raw_data/SRR_Acc_List.txt | while read ID;do 
    mkdir -p $HOME/gene_expression/trim/${ID} 
    trimmomatic PE -threads 6 -phred33 $HOME/gene_expression/raw_data/${ID}/${ID}_1.fastq $HOME/gene_expression/raw_data/${ID}/${ID}_2.fastq $HOME/gene_expression/trim/${ID}/${ID}_1.fq.gz $HOME/gene_expression/trim/${ID}/${ID}_1_unpaired.fq.gz $HOME/gene_expression/trim/${ID}/${ID}_2.fq.gz $HOME/gene_expression/trim/${ID}/${ID}_2_unpaired.fq.gz LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:30 ILLUMINACLIP:~/software/miniconda3/envs/placental_mosaicism/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa:2:30:10:1:true
done


## reads were aligned to the the hs1 NCBI gtf using STAR
mkdir $HOME/CHM13/star_index
STAR --runThreadN 6 --runMode genomeGenerate --genomeDir $HOME/CHM13/star_index_150 --genomeFastaFiles $HOME/CHM13/chm13v2.0.fa --sjdbGTFfile $HOME/CHM13/hs1.ncbiRefSeq.gtf --sjdbOverhang 150

cat $HOME/gene_expression/raw_data/SRR_Acc_List.txt | while read ID;do 
    mkdir -p $HOME/gene_expression/star/${ID}
    STAR --twopassMode Basic --genomeDir $HOME/CHM13/star_index_150 --runThreadN 6 --readFilesIn $HOME/gene_expression/trim/${ID}/${ID}_1.fq.gz $HOME/gene_expression/trim/${ID}/${ID}_2.fq.gz --readFilesCommand zcat --outFileNamePrefix $HOME/gene_expression/star/${ID}/${ID} --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 6 --quantMode TranscriptomeSAM GeneCounts --sjdbGTFfile $HOME/CHM13/hs1.ncbiRefSeq.gtf --outReadsUnmapped Fastx
done


## a transcriptome assembly was reconstructed using StringTie
cat $HOME/gene_expression/raw_data/SRR_Acc_List.txt | grep -v SRR24900365 | while read ID;do
    stringtie -p 8 -G $HOME/CHM13/hs1.ncbiRefSeq.gtf -o $HOME/gene_expression/stringtie/${ID}.gtf -i $HOME/gene_expression/star/${ID}/${ID}Aligned.sortedByCoord.out.bam
done

ls $HOME/gene_expression/stringtie/SRR*.gtf > $HOME/gene_expression/stringtie/merge_list.txt
stringtie --merge -p 8 -G $HOME/CHM13/hs1.ncbiRefSeq.gtf -o $HOME/gene_expression/stringtie/stringtie_merged.gtf $HOME/gene_expression/stringtie/merge_list.txt


## transcript tpm (transcripts per million) matrix was computed and spuriously expressed transcripts were removed (TPM <1 in at least half of samples)
mkdir $HOME/gene_expression/ballgown
cat $HOME/gene_expression/raw_data/SRR_Acc_List.txt | grep -v SRR24900365 | while read ID;do
    stringtie -e -B -p 8 -G $HOME/gene_expression/stringtie/stringtie_merged.gtf -o $HOME/gene_expression/ballgown/${ID}/${ID}.gtf $HOME/gene_expression/star/${ID}/${ID}Aligned.sortedByCoord.out.bam
done

prepDE_tpm.py -i $HOME/gene_expression/ballgown -l 150

cat $HOME/gene_expression/ballgown/transcript_tpm_matrix.csv |tr -t "," "\t" > $HOME/gene_expression/ballgown/transcript_tpm_matrix.txt

awk -F'\t' 'NR==1 {print; next} {count=0; for(i=2;i<=NF;i++) if($i > 1) count++; if(count >= NF/2) print}' $HOME/gene_expression/ballgown/transcript_tpm_matrix.txt|awk 'NR==1 {print; next} {sum=0; for(i=2;i<=NF;i++) sum+=$i; print $1"\t"sum/8}'|awk '$2>1' > $HOME/gene_expression/ballgown/transcript_tpm_matrix_all.txt


## the transcript with highest expression level of each gene was retained to generate a placenta-specific GTF.
awk -F "\t" '$9 ~ /^gene_id/ {print $9}' $HOME/gene_expression/stringtie/stringtie_merged.gtf|awk -F '[ ;]' '{print $2"\t"$5}'|uniq|sed 's/"//g'|sort -V|uniq > $HOME/gene_expression/ballgown/gene_transcript.txt

cat $HOME/gene_expression/ballgown/gene_transcript.txt|cut -f 1|sort -V|uniq|while read gene;do awk '{if(FNR==NR){s[$1]++}else{if(s[$1]==1){print $0}}}' <(awk -v gene=$gene '$1==gene' $HOME/gene_expression/ballgown/gene_transcript.txt|cut -f 2) $HOME/gene_expression/ballgown/transcript_tpm_matrix_all.txt|awk '($2>max||max==""){max=$2;line=$1}END{print line}';done|sed '/^$/d' > $HOME/gene_expression/ballgown/transcript_highest_expression_tpm.all.txt

grep -wFf <(sort $HOME/gene_expression/ballgown/transcript_highest_expression_tpm.all.txt|uniq|awk '{print "transcript_id \""$0"\""}') $HOME/gene_expression/stringtie/stringtie_merged.gtf|sort -k1,1 -k4,4n -k5,5n > $HOME/gene_expression/placenta.gtf
