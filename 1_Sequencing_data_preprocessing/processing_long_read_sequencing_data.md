configfile: "processing_long_read_sequencing_data.yaml"

rule all:
    input:
        expand("$HOME/nanopore_processing/nanostat/{sample}/{sample}.txt",sample=config['samples']),
        expand("$HOME/nanopore_processing/nanostat/{sample}-nanofilt/{sample}.txt",sample=config['samples']),
        expand("$HOME/nanopore_processing/nanoQC/{sample}/nanoQC.html",sample=config['samples']),
        expand("$HOME/nanopore_processing/nanoQC/{sample}-nanofilt/nanoQC.html",sample=config['samples']),
        expand("$HOME/nanopore_processing/nanoplot/{sample}/NanoStats.txt",sample=config['samples']),
        expand("$HOME/nanopore_processing/nanoplot/{sample}/HistogramReadlength.pdf",sample=config['samples']),
        expand("$HOME/nanopore_processing/nanoplot/{sample}/LengthvsQualityScatterPlot_hex.pdf",sample=config['samples']),
        expand("$HOME/nanopore_processing/nanoplot/{sample}-nanofilt/NanoStats.txt",sample=config['samples']),
        expand("$HOME/nanopore_processing/nanoplot/{sample}-nanofilt/HistogramReadlength.pdf",sample=config['samples']),
        expand("$HOME/nanopore_processing/nanoplot/{sample}-nanofilt/LengthvsQualityScatterPlot_hex.pdf",sample=config['samples']),
        expand("$HOME/nanopore_processing/nanofilt/{sample}.nanofilt.fq.gz",sample=config['samples']),
        expand("$HOME/nanopore_processing/minimap2/{sample}.bam",sample=config['samples']),
        expand("$HOME/nanopore_processing/minimap2/{sample}.bam.bai",sample=config['samples'])
        

rule nanostat:
    input:
        lambda wildcards: config['samples'][wildcards.sample][0]
    output:
        "$HOME/nanopore_processing/nanostat/{sample}/{sample}.txt"
    shell:
        "NanoStat --fastq {input} -t 8 -n {output} -o $HOME/nanopore_processing/nanostat/{wildcards.sample}"

rule nanoQC:
    input:
        lambda wildcards: config['samples'][wildcards.sample][0]
    output:
        "$HOME/nanopore_processing/nanoQC/{sample}/nanoQC.html"
    shell:
        "nanoQC {input} -o $HOME/nanopore_processing/nanoQC/{wildcards.sample}"
        
        
rule nanoplot:
    input:
        lambda wildcards: config['samples'][wildcards.sample][0]
    output:
        "$HOME/nanopore_processing/nanoplot/{sample}/NanoStats.txt",
        "$HOME/nanopore_processing/nanoplot/{sample}/HistogramReadlength.pdf",
        "$HOME/nanopore_processing/nanoplot/{sample}/LengthvsQualityScatterPlot_hex.pdf"
    shell:
        "NanoPlot --outdir $HOME/nanopore_processing/nanoplot/{wildcards.sample} --fastq {input} --loglength --format pdf --plots hex dot"
        

rule nanofilt:
    input:
        lambda wildcards: config['samples'][wildcards.sample][0]
    output:
        "$HOME/nanopore_processing/nanofilt/{sample}.nanofilt.fq.gz"
    shell:
        "gunzip -c {input} | NanoFilt -q 7 -l 1000 --headcrop 50 --tailcrop 20| gzip > $HOME/nanopore_processing/nanofilt/{wildcards.sample}.nanofilt.fq.gz"


rule nanofilt_stat:
    input:
        "$HOME/nanopore_processing/nanofilt/{sample}.nanofilt.fq.gz"
    output:
        "$HOME/nanopore_processing/nanostat/{sample}-nanofilt/{sample}.txt"
    shell:
        "NanoStat --fastq {input} -n {output} -o $HOME/nanopore_processing/nanostat/{wildcards.sample}-nanofilt"
        
        
rule nanofilt_QC:
    input:
        "$HOME/nanopore_processing/nanofilt/{sample}.nanofilt.fq.gz"
    output:
        "$HOME/nanopore_processing/nanoQC/{sample}-nanofilt/nanoQC.html"
    shell:
        "nanoQC {input} -o $HOME/nanopore_processing/nanoQC/{wildcards.sample}-nanofilt"


rule nanoplotClean:
    input:
        "$HOME/nanopore_processing/nanofilt/{sample}.nanofilt.fq.gz"
    output:
        "$HOME/nanopore_processing/nanoplot/{sample}-nanofilt/NanoStats.txt",
        "$HOME/nanopore_processing/nanoplot/{sample}-nanofilt/HistogramReadlength.pdf",
        "$HOME/nanopore_processing/nanoplot/{sample}-nanofilt/LengthvsQualityScatterPlot_hex.pdf"
    shell:
        "NanoPlot --outdir $HOME/nanopore_processing/nanoplot/{wildcards.sample}-nanofilt --fastq {input} --loglength --format pdf --plots hex dot"
        
        
rule minimap2:
    input:
        "$HOME/nanopore_processing/nanofilt/{sample}.nanofilt.fq.gz"
    output:
        "$HOME/nanopore_processing/minimap2/{sample}.bam"
    shell:
        "minimap2 -ax map-ont --MD -t 8 $HOME/CHM13/chm13v2.0.fa {input} | samtools sort -m 2G -O bam -o {output} -"
        

rule index:
    input:
        "$HOME/nanopore_processing/minimap2/{sample}.bam"
    output:
        "$HOME/nanopore_processing/minimap2/{sample}.bam.bai"
    shell:
        "samtools index -@ 8 {input} > {output}"
