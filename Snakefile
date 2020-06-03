configfile: "config.yaml"

import pandas as pd

data =  pd.read_table("data_table.tsv").set_index(['samples'], drop=False)

SAMPLE =  ['GSM787552_CD14_H3K36me3', 'GSM787528_CD14_H3K36me3', 'GSM787529_CD14_H3K36me3', 'GSM787524_CD14_H3K36me3','GSM787515_CD14_H3K36me3', 'GSM787513_CD14_H3K36me3']


rule all:
    input:
        expand("qc/fastqc/{samples}.html", samples = SAMPLE),
        expand("qc/fastqc/{samples}_fastqc.zip", samples = SAMPLE),
        "qc/multiqc/reads.html",
        expand("indexes/{genome}/{organism}.1.bt2", genome = config['genome'], organism = config['organism']),
        expand("indexes/{genome}/{organism}.2.bt2", genome = config['genome'], organism = config['organism']),
        expand("indexes/{genome}/{organism}.3.bt2", genome = config['genome'], organism = config['organism']),
        expand("indexes/{genome}/{organism}.4.bt2", genome = config['genome'], organism = config['organism']),
        expand("indexes/{genome}/{organism}.rev.1.bt2", genome = config['genome'], organism = config['organism']),
        expand("indexes/{genome}/{organism}.rev.2.bt2", genome = config['genome'], organism = config['organism']),
        expand("mapped/{samples}.sam", samples = SAMPLE),
        expand("mapped/{samples}.bam", samples = SAMPLE),
        expand("mapped/{samples}_sorted.bam", samples = SAMPLE),
        expand("mapped/{samples}_sorted.bam.bai", samples = SAMPLE),
        expand("bw/{samples}.bw", samples = SAMPLE),
        "qc/multiqc/bams.html",
        "result.tar.gz"
       
rule fastqc:
    input:
        lambda wildcards: data['File'][wildcards.samples]
    output:
        html="qc/fastqc/{samples}.html",
        zip="qc/fastqc/{samples}_fastqc.zip" 
    params: ""
    log:
        "logs/fastqc/{samples}.log"
    wrapper:
        "0.59.2/bio/fastqc"
rule multiqc:
    input:
        expand("qc/fastqc/{samples}_fastqc.zip", samples = SAMPLE)
    output:
        "qc/multiqc/reads.html"
    shell:
        "multiqc -d {input} -n {output}"

rule wget:
    output: expand("indexes/{genome}/{genome}.fa.gz", genome = config['genome'])
    priority: 1
    params:
        link = config['genome']
    shell: "wget -O {output} http://hgdownload.soe.ucsc.edu/goldenPath/{params.link}/bigZips/{params.link}.fa.gz" 

rule bowtie2Build:
    input: expand("indexes/{genome}/{genome}.fa.gz", genome = config['genome'])
    params:
        basename = expand("indexes/{genome}/{organism}", genome = config['genome'], organism = config['organism'])
    output:
        expand("indexes/{genome}/{organism}.1.bt2", genome = config['genome'], organism = config['organism']),
        expand("indexes/{genome}/{organism}.2.bt2", genome = config['genome'], organism = config['organism']),
        expand("indexes/{genome}/{organism}.3.bt2", genome = config['genome'], organism = config['organism']),
        expand("indexes/{genome}/{organism}.4.bt2", genome = config['genome'], organism = config['organism']),
        expand("indexes/{genome}/{organism}.rev.1.bt2", genome = config['genome'], organism = config['organism']),
        expand("indexes/{genome}/{organism}.rev.2.bt2", genome = config['genome'], organism = config['organism'])
    priority: 1
    conda:
        "envs/bowtie2.yaml"
    shell: "bowtie2-build {input} {params.basename}"

rule bowtie_aln:
    input: 
        lambda wildcards: data['File'][wildcards.samples]
    params:
        basename=expand("indexes/{genome}/{organism}", genome = config['genome'], organism = config['organism'])
    output: "mapped/{samples}.sam"
    log: "logs/bowtie2/{samples}.log"
    conda:
        "envs/bowtie2.yaml"
    shell:
        "bowtie2 -x {params.basename} -U {input} 2> {log} > {output}"

rule sam2bam:
    input: "mapped/{samples}.sam"
    output:
        bam = "mapped/{samples}.bam", 
        bam_sort = "mapped/{samples}_sorted.bam",
        bai = "mapped/{samples}_sorted.bam.bai"
    run:
        shell("samtools view -Sb {input} > {output.bam}"),
        shell("samtools sort {output.bam} -o {output.bam_sort}"),
        shell("samtools index {output.bam_sort}")

rule bigwig:
    input:
        "mapped/{samples}_sorted.bam"
    output:
        "bw/{samples}.bw"
    conda:
        "envs/bigwig.yaml"
    shell: "bamCoverage -b {input} -o {output}"

rule multiqc_bam:
    input:
        expand("logs/bowtie2/{samples}.log", samples = SAMPLE)
    output:
        "qc/multiqc/bams.html"
    shell: "multiqc -d {input} -n {output}"

rule tar:
        output: "result.tar.gz"
        shell: "tar -zcvf result.tar.gz bw qc/multiqc"
