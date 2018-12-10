import yaml
import re

##========= Globals ==========
configfile: 'config.yaml'

## Set samples
FILES = yaml.load(open(config['sampleYaml']))
SAMPLES = sorted(FILES.keys())

## set output Dir 
WORKDIR = config["workDir"]

## Set reference file
DNA = config["dna"]
GTF = config["gtf"]
cDNA = config["cDNA"]
HISAT2_INDEX_PREFIX = config["genome_hisat2_index"]

##======== Rules ============
rule all:
    input:
        expand( WORKDIR + "Step0.Prepare/{sample}.R1.fq.gz", sample=SAMPLES ),
        expand( WORKDIR + "Step0.Prepare/{sample}.R2.fq.gz", sample=SAMPLES ),
        expand( WORKDIR + "Step1.fastqFilter/{sample}/{sample}.R1.fq.gz", sample=SAMPLES ),
        expand( WORKDIR + "Step1.fastqFilter/{sample}/{sample}.R2.fq.gz", sample=SAMPLES ),
        expand( WORKDIR + "Step1.fastqFilter/{sample}/{sample}.json", sample=SAMPLES ),
        expand( WORKDIR + "Step1.fastqFilter/{sample}/{sample}.html", sample=SAMPLES ),
        expand( WORKDIR + "Step2.hisat2Align/{sample}.bam", sample=SAMPLES),
        expand( WORKDIR + "Step2.sorted/{sample}.sorted.bam", sample=SAMPLES),
        expand( WORKDIR + "Step3.stringtieAssemble/{sample}/{sample}.gtf", sample=SAMPLES ),
        WORKDIR + "Step4.stringtieMerge/merged_list.txt",
        WORKDIR + "Step4.stringtieMerge/stringtie_merged.gtf",
        expand( WORKDIR + "Step5.stringtieQuant/{sample}/{sample}.gtf", sample=SAMPLES ),
        WORKDIR + "Step6.StringTieDE/sample_gtf.list"
 

## ======== Step 0: Prepare Work ========
rule ReName:
    input:
        lambda wildcards:FILES[wildcards.sample]
    output:
       R1 =  WORKDIR + "Step0.Prepare/{sample}.R1.fq.gz",
       R2 =  WORKDIR + "Step0.Prepare/{sample}.R2.fq.gz"
    run:
        shell("ln -sf {input[0]} {output.R1}")
        shell("ln -sf {input[1]} {output.R2}")


##======== Step1 raw fastq filter ========
rule FastqFilter:
    input:
        R1 =  WORKDIR + "Step0.Prepare/{sample}.R1.fq.gz",
        R2 =  WORKDIR + "Step0.Prepare/{sample}.R2.fq.gz"
    output:
        R1 = WORKDIR + "Step1.fastqFilter/{sample}/{sample}.R1.fq.gz",
        R2 = WORKDIR + "Step1.fastqFilter/{sample}/{sample}.R2.fq.gz",
        json = WORKDIR + "Step1.fastqFilter/{sample}/{sample}.json",
        html = WORKDIR + "Step1.fastqFilter/{sample}/{sample}.html"
    log:
        WORKDIR + "logs/Step1.fastqFilter/{sample}.fastqFilter.log"
    benchmark:
        WORKDIR + "benchmark/Step1.fastqFilter/{sample}.benchmark"
    params:
        " --detect_adapter_for_pe --qualified_quality_phred 20 --unqualified_percent_limit 5 "
        " --n_base_limit 5 --length_required 50 --correction "
    threads:
        8
    run:
        shell("fastp {params} -w {threads} -i {input.R1} -I {input.R2} -o {output.R1} "
        "-O {output.R2} -j {output.json} -h {output.html} 2> {log}")


## ======== Step2 hisat alignment ========
rule hisat2Align:
    input:
        R1 = WORKDIR + "Step1.fastqFilter/{sample}/{sample}.R1.fq.gz",
        R2 = WORKDIR + "Step1.fastqFilter/{sample}/{sample}.R2.fq.gz"
    output:
        temp(WORKDIR + "Step2.hisat2Align/{sample}.bam")
    log:
        WORKDIR + "logs/Step2.hisat2Align/{sample}.align.logs"
    benchmark:
        WORKDIR + "benchmark/Step2.hisat2Align/{sample}.benchmark"
    threads:
        8
    params:
        "--dta"
    shell:
         "hisat2 {params} -p {threads} -x {HISAT2_INDEX_PREFIX} -1 {input.R1} -2 {input.R2} 2>{log} |samtools view -Sb -@ {threads} -m 5G -o {output}"

## ======== Step 2.2 sort bam ========
rule sortBam:
    input:
        bam = WORKDIR + "Step2.hisat2Align/{sample}.bam"
    output:
        protected(WORKDIR + "Step2.sorted/{sample}.sorted.bam")
    threads:
        8
    resources: 
        mem_mb=5000
    shell:
        "samtools sort -@ {threads} -m 5G -o {output} {input.bam}"

## ======== Step3 stringtie assemble ========
rule stringtieAssemble:
    input:
        bam = WORKDIR + "Step2.sorted/{sample}.sorted.bam"
    output: 
        WORKDIR + "Step3.stringtieAssemble/{sample}/{sample}.gtf"
    log:
        WORKDIR + "logs/Step3.stringtieAssemble/{sample}.stringtieAssemble.log"
    benchmark:
        WORKDIR + "benchmark/Step3.stringtieAssembl{sample}.benchmark"
    params:
        "-m 250 -j 3 -c 2.5 -B"
    threads:
        8
    run:
       shell( "stringtie {input.bam} {params} -p {threads} -G {GTF} -o {output} 2>{log}" )


## ======= Step4.1 make StringTie merge list ======
rule stringtieMergeList:
     input: 
         expand(WORKDIR + "Step3.stringtieAssemble/{sample}/{sample}.gtf", sample = SAMPLES)
     output: 
         list = WORKDIR + "Step4.stringtieMerge/merged_list.txt"
     run:
        with open(output.list, 'w') as f:
            for gtf in input:
                print(gtf, file=f)

## ======== Step4 stringtie merge ========
rule stringtieMerge:
    input:
        list = WORKDIR + "Step4.stringtieMerge/merged_list.txt"
    output:
        WORKDIR + "Step4.stringtieMerge/stringtie_merged.gtf"
    params:
        "-m 50 -c 3 "
    threads:
        8
    run:
       shell( "stringtie --merge {params} -p {threads} -G {GTF} -o {output} {input.list}")


## ======= Step 5 stringtie quant =======
rule stringtieQuant:
    input:
        merged_gtf = WORKDIR + "Step4.stringtieMerge/stringtie_merged.gtf",
        bam = WORKDIR + "Step2.sorted/{sample}.sorted.bam"
    output:
        gtf = WORKDIR + "Step5.stringtieQuant/{sample}/{sample}.gtf"
    log:
        WORKDIR + "logs/Step5.stringtieQuant/{sample}.log" 
    benchmark:
        WORKDIR + "benchmark/Step5.stringtieQuant/{sample}.benchmark"
    params:
        "-e -B"
    threads: 
        8
    run:
        shell("stringtie {params} -p {threads} -G {input.merged_gtf} -o {output.gtf} {input.bam} 2> {log}")

## ======== Step 6.1 Difference Expression Analysis with prepDE.py ========
rule produce_sample_express_list:
    input:
        expand(WORKDIR + "Step5.stringtieQuant/{sample}/{sample}.gtf", sample=SAMPLES)
    output:
        list = WORKDIR + "Step6.StringTieDE/sample_gtf.list"
    run:
        with open(output.list, 'w') as f:
            for gtf in input:
                sample = re.search("stringtieQuant\/(\S+)\/(\S+).gtf", gtf)
                out = sample.group(1) + "\t" + gtf
                print(out, file=f)

## ======== Step 6.2 Fetch Gene Expression Matrix ========
rule prepDE:
    input:
        list = WORKDIR + "Step6.StringTieDE/sample_gtf.list"
    output:
        gene = WORKDIR + "Step6.StringTieDE/gene_count_matrix.csv",
        transcript = WORKDIR + "Step6.StringTieDE/transcript_count_matrix.csv"
    run:
        shell("prepDE.py -i {input.list} -g {output.gene} -t {output.transcript}")
