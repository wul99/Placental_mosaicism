configfile: "processing_short_read_sequencing_data.yaml"

rule all:
    input:
        expand("$HOME/raw_data_cutadapter/{sample}/{sample}_{num}_clean.fq.gz",sample=config['samples'],num=['1', '2']),
        expand("$HOME/raw_data_cutadapter_bwa/{sample}.sorted.bam",sample=config['samples']),
        expand("$HOME/raw_data_cutadapter_bwa_rmdup/{sample}.rmdup.sorted.bam",sample=config['samples']),
        expand("$HOME/raw_data_cutadapter_bwa_rmdup/{sample}.rmdup.reID.sorted.bam",sample=config['samples']),
        expand("$HOME/raw_data_cutadapter_bwa_rmdup_bqsr_byCHM13/{sample}.rmdup.bqsr.recal_data.table",sample=config['samples']),
        expand("$HOME/raw_data_cutadapter_bwa_rmdup_bqsr_byCHM13/{sample}.rmdup.bqsr.bam",sample=config['samples']),
        expand("$HOME/raw_data_cutadapter_bwa_quality_byCHM13/{sample}.tab",sample=config['samples']),
        expand("$HOME/raw_data_cutadapter_bwa_quality_byCHM13/{sample}.rmdup.txt",sample=config['samples']),
        expand("$HOME/raw_data_cutadapter_bwa_quality_byCHM13/{sample}.rmdup.bqsr.txt",sample=config['samples']),
        expand("$HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/{sample}.bam",sample=config['samples']),
        expand("$HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/{sample}.bam.bai",sample=config['samples']),
        expand("$HOME/raw_data_cutadapter_bwa_quality_byCHM13/{sample}.rmdup.bqsr.filter_Q30.txt",sample=config['samples'])

rule cutadapter:
    input:
        lambda wildcards: config['samples'][wildcards.sample][0],
        lambda wildcards: config['samples'][wildcards.sample][1]
    output:
        "$HOME/raw_data_cutadapter/{sample}/{sample}_1_clean.fq.gz",
        "$HOME/raw_data_cutadapter/{sample}/{sample}_2_clean.fq.gz"
    shell:
        "SOAPnuke filter -n 0.1 -m 30 -l 12 -T 8 -J -1 {input[0]} -2 {input[1]} -f AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA -r AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG -o $HOME/raw_data_cutadapter/{wildcards.sample} -C {wildcards.sample}_1_clean.fq.gz -D {wildcards.sample}_2_clean.fq.gz"

rule bwa:
    input:
        "$HOME/raw_data_cutadapter/{sample}/{sample}_1_clean.fq.gz",
        "$HOME/raw_data_cutadapter/{sample}/{sample}_2_clean.fq.gz"
    output:
        "$HOME/raw_data_cutadapter_bwa/{sample}.sorted.bam"
    shell:
        "~/software/miniconda3/envs/placental_mosaicism/bin/bwa-mem2 mem -t 8 -M $HOME/CHM13/chm13v2.0.fa {input[0]} {input[1]} | samtools view -bS | samtools sort -o {output}"

rule rm_dup:
    input:
        "$HOME/raw_data_cutadapter_bwa/{sample}.sorted.bam"
    output:
        "$HOME/raw_data_cutadapter_bwa_rmdup/{sample}.rmdup.sorted.bam"
    shell:
        "sh $HOME/rm_duplicates.sh {input} {output} {wildcards.sample}"

rule reID:
    input:
        "$HOME/raw_data_cutadapter_bwa_rmdup/{sample}.rmdup.sorted.bam",
        "$HOME/dataname_all.txt"
    output:
        "$HOME/raw_data_cutadapter_bwa_rmdup/{sample}.rmdup.reID.sorted.bam"
    shell:
        "picard AddOrReplaceReadGroups I={input[0]} O={output} RGID={wildcards.sample} RGLB=$(grep -w {wildcards.sample} {input[1]}|cut -f 3) RGPL=DNBSEQ RGPU=$(grep -w {wildcards.sample} {input[1]}|cut -f 4) RGSM=$(grep -w {wildcards.sample} {input[1]}|cut -f 1) VALIDATION_STRINGENCY=LENIENT"

rule BaseRecalibrator:
    input:
        "$HOME/raw_data_cutadapter_bwa_rmdup/{sample}.rmdup.reID.sorted.bam"
    output:
        "$HOME/raw_data_cutadapter_bwa_rmdup_bqsr_byCHM13/{sample}.rmdup.bqsr.recal_data.table"
    shell:
        "gatk BaseRecalibrator --tmp-dir $HOME/raw_data_cutadapter_bwa_rmdup_bqsr_byCHM13/tmp -R $HOME/CHM13/chm13v2.0.fa -I {input} --known-sites $HOME/CHM13/1000G_phase1.indels.hs1.sites.sorted.vcf --known-sites $HOME/CHM13/Mills_and_1000G_gold_standard.indels.hs1.sites.sorted.vcf --known-sites $HOME/CHM13/dbsnp_138.hs1.filter.sorted.vcf -O {output}"

rule ApplyBQSR:
    input:
        "$HOME/raw_data_cutadapter_bwa_rmdup/{sample}.rmdup.reID.sorted.bam",
        "$HOME/raw_data_cutadapter_bwa_rmdup_bqsr_byCHM13/{sample}.rmdup.bqsr.recal_data.table"
    output:
        "$HOME/raw_data_cutadapter_bwa_rmdup_bqsr_byCHM13/{sample}.rmdup.bqsr.bam"
    shell:
        "gatk ApplyBQSR -R $HOME/CHM13/chm13v2.0.fa -I {input[0]} -O {output} -bqsr {input[1]} --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 --add-output-sam-program-record"

rule samtools_filter:
    input:
        "$HOME/raw_data_cutadapter_bwa_rmdup_bqsr_byCHM13/{sample}.rmdup.bqsr.bam"
    output:
        "$HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/{sample}.bam"
    shell:
        "samtools view -@ 8 -q 30 -F 3340 {input} -b -o {output}"

rule samtools_index:
    input:
        "$HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/{sample}.bam"
    output:
        "$HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/{sample}.bam.bai"
    shell:
        "samtools index -@ 8 {input} {output}"

rule flagstat:
    input:
        "$HOME/raw_data_cutadapter_bwa/{sample}.sorted.bam"
    output:
        "$HOME/raw_data_cutadapter_bwa_quality_byCHM13/{sample}.tab"
    shell:
        "samtools flagstat -@ 8 {input} > {output}"

rule rmdup_flagstat:
    input:
        "$HOME/raw_data_cutadapter_bwa_rmdup/{sample}.rmdup.sorted.bam"
    output:
        "$HOME/raw_data_cutadapter_bwa_quality_byCHM13/{sample}.rmdup.txt"
    shell:
        "samtools flagstat -@ 8 {input} > {output}"

rule bqsr_flagstat:
    input:
        "$HOME/raw_data_cutadapter_bwa_rmdup_bqsr_byCHM13/{sample}.rmdup.bqsr.bam"
    output:
        "$HOME/raw_data_cutadapter_bwa_quality_byCHM13/{sample}.rmdup.bqsr.txt"
    shell:
        "samtools flagstat -@ 8 {input} > {output}"

rule filter_flagstat:
    input:
        "$HOME/raw_data_cutadapter_bwa_rmdup_bqsr_Q30_byCHM13/{sample}.bam"
    output:
        "$HOME/raw_data_cutadapter_bwa_quality_byCHM13/{sample}.rmdup.bqsr.filter_Q30.txt"
    shell:
        "samtools flagstat -@ 8 {input} > {output}"
