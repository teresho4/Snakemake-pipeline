
SAMPLES = ["SRR787528.chr15", "SRR787552.chr15"]

rule all:
    input:
        "indexes/human.1.bt2",
        "indexes/human.2.bt2",
        "indexes/human.3.bt2",
        "indexes/human.4.bt2",
        "indexes/human.rev.1.bt2",
        "indexes/human.rev.2.bt2"

rule fastqc:
    input:
        "reads/{sample}.fastq"
    output:
        "qc/fastqc/{sample}_fastqc.html",
        "qc/fastqc/{sample}_fastqc.zip"
    log:
        "logs/fastqc/{sample}.log"
    shell:
        "fastqc {input} -q -o qc/fastqc/"

rule multiqc:
    input:
        expand("qc/fastqc/{sample}_fastqc.html", sample = SAMPLES)
    output:
        "qc/multiqc_report.html"
    wrapper:
        "0.31.1/bio/multiqc"

rule wget:
    output:
 	"chr15/chr15.fa.gz"
    shell:
        "wget -O ./chr15/chr15.fa.gz 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr15.fa.gz'"

rule bowtie2Build:
    input: "chr15/chr15.fa.gz"
    params:
        basename = "./indexes/human"
    output:
        output1="indexes/human.1.bt2",
        output2="indexes/human.2.bt2",
        output3="indexes/human.3.bt2",
        output4="indexes/human.4.bt2",
        outputrev1="indexes/human.rev.1.bt2",
        outputrev2="indexes/human.rev.2.bt2"
    log:
	"logs/bowtie2/bowtie-build.log"
    shell: "bowtie2-build {input} {params.basename}"



