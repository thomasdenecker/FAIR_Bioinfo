"""
Author: T. Denecker & C. Toffano-Nioche
Affiliation: I2BC
Aim: A workflow to process RNA-Seq.
Organism : O. tauri

Date: Jan 2019
Run: snakemake -s Snakefile
Processus:
    - Quality control : fastqc
    - Indexing of the reference genome : bowtie2 BUILD
    - Alignment of reads with the reference genome
    - Manipulation of alignment files
        - binary conversion,
        - sorting
        - indexing
    - Counting the number of reads aligned on each gene

Latest modification:
  - todo
"""

##----------------------------------------------------------------------------##
## Working directory                                                          ##
##----------------------------------------------------------------------------##

BASE_DIR = "/home/rstudio"
WDIR = BASE_DIR + "/Project"
workdir: WDIR

##----------------------------------------------------------------------------##
## Variables declaration
## Declaring some variables used by topHat and other tools...
## (GTF file, INDEX, chromosome length)
##----------------------------------------------------------------------------##

GFF   = WDIR + "/annotations/O.tauri_annotation.gff"
GENOME = WDIR + "/genome/O.tauri_genome.fna"

##----------------------------------------------------------------------------##
## The list of samples to be processed
##----------------------------------------------------------------------------##
SAMPLES, = glob_wildcards("./samples/{smp}.fastq.gz")
NB_SAMPLES = len(SAMPLES)

##----------------------------------------------------------------------------##
## Rule final
##----------------------------------------------------------------------------##

rule final:
    input:
        expand("fastqc/{smp}/{smp}_fastqc.zip", smp=SAMPLES),
        expand("htseq/count_{smp}.txt", smp=SAMPLES),
        expand("graphics/graphics-{smp}.pdf", smp=SAMPLES),
        expand("countTable.txt")

##----------------------------------------------------------------------------##
## Fastqc
##----------------------------------------------------------------------------##

rule fastqc:
        input:  "samples/{smp}.fastq.gz"
        output: "fastqc/{smp}/{smp}_fastqc.zip"
        message: """--- Quality check of raw data with Fastqc."""
        shell: """
        fastqc {input} --outdir  fastqc/{wildcards.smp}
        """
##----------------------------------------------------------------------------##
## BOWTIE2 BUILD
##----------------------------------------------------------------------------##

rule bowtie2Build:
    input: GENOME
    params:
        basename= "reference/reference"
    output:
        output1="reference/reference.1.bt2",
        output2="reference/reference.2.bt2",
        output3="reference/reference.3.bt2",
        output4="reference/reference.4.bt2",
        outputrev1="reference/reference.rev.1.bt2",
        outputrev2="reference/reference.rev.2.bt2"
    message: """--- Indexation of the reference genome."""
    shell: "bowtie2-build {input} {params.basename}"

##----------------------------------------------------------------------------##
## BOWTIE2
##----------------------------------------------------------------------------##

rule bowtie2:
    input:
        file = "samples/{smp}.fastq.gz",
         bt2 = "reference/reference.rev.2.bt2"
    params:
        index = 'reference/reference'
    output:
        sam = "bowtie2/bowtie-{smp}.sam",
        out = "bowtie2/bowtie-{smp}.out"
    message: """--- Alignment of reads with the reference genome."""
    shell: 'bowtie2 -x {params.index} -U {input.file} -S {output.sam} > {output.out}'

##----------------------------------------------------------------------------##
## Samtools - View
##----------------------------------------------------------------------------##
rule samtoolsView:
    input:
        "bowtie2/bowtie-{smp}.sam"
    output:
        "samtools/bowtie-{smp}.bam"
    message: """--- Binary conversion of aligned reads."""
    shell: """
    samtools view -b {input} > {output}
    """

##----------------------------------------------------------------------------##
## Samtools - Sort and index
##----------------------------------------------------------------------------##
rule samtoolsSortIndex:
    input:
        "samtools/bowtie-{smp}.bam"
    output:
        "samtools/bowtie-{smp}.sorted.bam"
    message: """--- Sorting and Indexingof aligned reads."""
    shell: """
    samtools sort {input} -o {output}
    samtools index {output}
    """
##----------------------------------------------------------------------------##
## Hteseq-count
##----------------------------------------------------------------------------##

rule htseqCount:
    input:
        "samtools/bowtie-{smp}.sorted.bam"
    output:
        "htseq/count_{smp}.txt"
    params:
        GFF
    message: """--- Count."""
    shell: """
    htseq-count --stranded=no --type='gene' --idattr='ID' --order=name --format=bam {input} {params} > {output}
    """

##----------------------------------------------------------------------------##
## Rscript
##----------------------------------------------------------------------------##

rule graphics:
    input:
        "htseq/count_{smp}.txt"
    output:
        "graphics/graphics-{smp}.pdf"
    params:
        GFF
    message: "--- Histogram"
    shell: """
    Rscript ../R-code/graphics.R {input}
    """

##----------------------------------------------------------------------------##
## Rscript
##----------------------------------------------------------------------------##

rule CountTable:
    input:
        "../conditions.txt"
    output:
        "countTable.txt"
    params:
        GFF
    message: "--- Create count table"
    shell: """
    Rscript ../R-code/countTable.R
    """
