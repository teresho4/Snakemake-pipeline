
configfile: "config.yaml"

rule all:
    input:
        "indexes/chr15/human.1.bt2",
        "indexes/chr15/human.2.bt2",
        "indexes/chr15/human.3.bt2",
        "indexes/chr15/human.4.bt2",
        "indexes/chr15/human.rev.1.bt2",
        "indexes/chr15/human.rev.2.bt2",
	expand("mapped/{sample}.sam", sample = config["samples"]),
	expand("mapped/{sample}.bam", sample = config["samples"]),
	expand("mapped/{sample}_sorted.bam", sample = config["samples"]),
	expand("mapped/{sample}_sorted.bam.bai", sample = config["samples"]),
	expand("qc/fastqc/{sample}_fastqc.html", sample = config["samples"]),
    	expand("qc/fastqc/{sample}_fastqc.zip", sample = config["samples"]),
	"qc/multiqc/reads.html",
	expand("bw/{sample}.bw", sample = config["samples"]),
        expand("logs/bowtie2/{sample}.log", sample = config["samples"]),
        "qc/multiqc/bams.html",
	"result.tar.gz"

rule fastqc:
    input:
        "reads/{sample}.fastq"
    output:
        "qc/fastqc/{sample}_fastqc.html",
        "qc/fastqc/{sample}_fastqc.zip"
    shell:
        "fastqc {input} -q -o qc/fastqc/"


rule multiqc:
    input:
        expand("qc/fastqc/{sample}_fastqc.zip", sample = config["samples"])
    output:
        "qc/multiqc/reads.html"
    shell:
        "multiqc -d {input} -n {output}"

rule wget:
    output: "indexes/chr15/chr15.fa.gz"
    priority: 1
    shell: "wget -O {output} ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr15.fa.gz"

rule bowtie2Build:
    input: "indexes/chr15/chr15.fa.gz"
    params:
        basename = "indexes/chr15/human"
    output:
        output1="indexes/chr15/human.1.bt2",
        output2="indexes/chr15/human.2.bt2",
        output3="indexes/chr15/human.3.bt2",
        output4="indexes/chr15/human.4.bt2",
        outputrev1="indexes/chr15/human.rev.1.bt2",
        outputrev2="indexes/chr15/human.rev.2.bt2"
    priority: 2
    shell: "bowtie2-build {input} {params.basename}"

rule bowtie_aln:
    input: "reads/{sample}.fastq"
    params:
        basename="indexes/chr15/human"
    output: "mapped/{sample}.sam"
    log: "logs/bowtie2/{sample}.log"
    shell:
        "bowtie2 -x {params.basename} -U {input} 2> {log} > {output}"

rule multiqc_bam:
    input:
        expand("logs/bowtie2/{sample}.log", sample = config["samples"])
    output:
        "qc/multiqc/bams.html"
    shell: "multiqc -d {input} -n {output}"

rule sam2bam:
    input: "mapped/{sample}.sam"
    output: "mapped/{sample}.bam"
    shell:
        "samtools view -Sb {input} > {output}"

rule index_bam:
        input: "mapped/{sample}.bam"
        output: "mapped/{sample}_sorted.bam"
        shell: "samtools sort {input} -o {output}"
rule index:
   	input: "mapped/{sample}_sorted.bam"
	priority: 2
        output: "mapped/{sample}_sorted.bam.bai"
        shell: "samtools index {input}"

rule bigwig:
    input:
        "mapped/{sample}_sorted.bam"
    output:
        "bw/{sample}.bw"
    priority: 1
    shell: "bamCoverage -b {input} -o {output}"

rule tar:
        output: "result.tar.gz"
        shell: "tar -zcvf result.tar.gz bw qc"
