# To snake or to flow: comparison of Snakemake and Nextflow
---
Authors: [Lloyd Low](https://researchers.adelaide.edu.au/profile/wai.low)

Role: Senior Bioinformatician

Affiliation: University of Adelaide

*14/12/2023*

## Table of Contents
* [Workflow management system](#workflow)
* [Prerequisites](#pre)
* [Dataset](#data)
* [Shell scripting](#shell)
* [Snakemake](#snakemake)
* [Nextflow](#nextflow)
* [Acknowledgements](#Acknowledgements)

<a name="workflow"></a>
## Workflow management system
Computational pipelines are integral in bioinformatics, particularly for high-throughput sequenced data. While traditional Shell scripting is a common method, alternatives like Snakemake and Nextflow have gained popularity. Notably, the [verkko assembler](https://github.com/marbl/verkko) for assembling complete human telomere-to-telomere genomes uses a Snakemake pipeline whereas many workflows customised to the use of Oxford Nanopore Technologies long reads use [Nextflow pipelines](https://labs.epi2me.io/wfindex/). This workshop addresses common queries such as the choice between Shell scripting and workflow management systems like Snakemake or Nextflow. By demonstrating variant calling with both Shell scripts and these systems, the workshop provides insights into factors to consider, including scalability. Pros and cons of each method underscore the idea that the optimal approach depends on specific task requirements and features.

Why you should consider learning a proper workflow manager?
* Reproducibility ([think container](https://www.nextflow.io/docs/latest/container.html))
* Efficiency (think parallel executions)
* Ease of resuming partially completed runs

<a name="pre"></a>
## Prerequisites

* Knowledge in shell scripting is assumed for this workshop

<a name="data"></a>
## Dataset

<details>
<summary>
<i> Raw fastq files and a reference genome </i>
</summary>
<p>

### Raw fastq files

Sample1.R1.fastq.gz

Sample1.R2.fastq.gz

Sample2.R1.fastq.gz

Sample2.R2.fastq.gz

Sample3.R1.fastq.gz

Sample3.R2.fastq.gz

### Reference genome (only chr 28 and chr 29)

ARS-UCD2.0_demo.fa

ARS-UCD2.0_demo.fa.amb

ARS-UCD2.0_demo.fa.ann

ARS-UCD2.0_demo.fa.bwt

ARS-UCD2.0_demo.fa.fai

ARS-UCD2.0_demo.fa.pac

ARS-UCD2.0_demo.fa.sa

</p>
</details>

<a name="shell"></a>
## Shell scripting

This shell script will process the raw fastq files, map to a reference genome and call single nucleotide polymorphisms (SNPs).

The steps are

* Raw reads fastqc
* Clean and trim reads
* Fastqc on trimmed reads
* Mapping to a ref with BWA
* Mark duplicates with Picard or SamBamba
* Samtools stats, idxstats, flagstats on the BAM files
* Call variants using bcfools (OPTIONAL: split by chromosome)
* Multiqc report

```bash
#!/bin/bash

#conda activate insyb2023 #ensure all tools such as fastqc, seqkit etc are available

#define variables
taskcpus=8
read1=Sample1.R1.fastq.gz
read2=Sample1.R2.fastq.gz
sampleID=Sample1
read1trimmed=Sample1.R1_val_1.fq.gz
read2trimmed=Sample1.R2_val_2.fq.gz
ref=ARS-UCD2.0_demo.fa
avail_mem=3
bam=Sample1.sorted.bam
dedupbam=Sample1.markDup.bam

#ensure ref and raw fastq are here
# ln -s ../raw_fastq/Sample1.R1.fastq.gz ../raw_fastq/Sample1.R2.fastq.gz .
# ln -s ../reference/* .

#fastqc_raw
fastqc -t ${taskcpus} ${read1} ${read2}

#seqkit_stats
seqkit stats -j ${taskcpus} ${read1} ${read2} > ${sampleID}.seqkit.stats

#trim_galore
trim_galore --quality 30 \
    --length 75 \
    --clip_R1 1 \
    --clip_R2 1 \
    --three_prime_clip_R1 1 \
    --three_prime_clip_R2 1 \
    --paired \
    --cores ${taskcpus} \
    ${read1} ${read2}

#fastqc_trimmed
fastqc -t ${taskcpus} ${read1trimmed} ${read2trimmed}

#bwa_mapping
bwa mem -M \
        -t ${taskcpus} \
        -R "@RG\tID:${sampleID}\tPL:${sampleID}\tLB:${sampleID}\tSM:${sampleID}" \
        $ref \
        ${read1trimmed} ${read2trimmed} \
        | samtools sort --threads ${taskcpus} -o ${sampleID}.sorted.bam

#picard_markDuplicates
java -jar -Xmx${avail_mem}g /home/lloyd/Software/picard/build/libs/picard.jar MarkDuplicates \
                I=${bam} \
                O=${sampleID}.markDup.bam \
                M=${sampleID}.markDup.metrics.txt \
                VALIDATION_STRINGENCY=LENIENT \
                REMOVE_DUPLICATES=false \
                OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 \
                CREATE_INDEX=false
samtools index ${sampleID}.markDup.bam

#bam_stats_pre_dedup
samtools index ${bam}
samtools stats ${bam} > ${sampleID}.prededup.stats
samtools flagstat ${bam} > ${sampleID}.prededup.flagstat
samtools idxstats ${bam} > ${sampleID}.prededup.idxstats

#bam_stats_post_dedup
samtools stats ${dedupbam} > ${sampleID}.postdedup.stats
samtools flagstat ${dedupbam} > ${sampleID}.postdedup.flagstat
samtools idxstats ${dedupbam} > ${sampleID}.postdedup.idxstats

#split_bam
samtools idxstats ${dedupbam} | cut -f 1 | grep -v '*' > ${sampleID}.chromosomes.txt
while IFS= read -r line; do
    samtools view -b ${dedupbam} ${line} > ${sampleID}.${line}.bam ;
    samtools index ${sampleID}.${line}.bam
done < ${sampleID}.chromosomes.txt

#call_genotype
while IFS= read -r line; do
    chromosome=$line
    finalbamExample=${sampleID}.${chromosome}.bam

    bcftools mpileup -f ${ref} ${finalbamExample} | \
    bcftools call --threads ${taskcpus} -m --output-type z | \
    bcftools filter --threads ${taskcpus} --output-type z \
    -s LowQual -e 'QUAL<20 || DP<100' > ${sampleID}.${chromosome}.vcf.gz
    bcftools index ${sampleID}.${chromosome}.vcf.gz
    bcftools stats ${sampleID}.${chromosome}.vcf.gz > ${sampleID}.${chromosome}.vcf.stats

done < ${sampleID}.chromosomes.txt
```  

<a name="snakemake"></a>
## Snakemake

![Snakemake](snakemake.png)

I recommend [Snakemake tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html) and this paper by [Koster and Rahmann 2012](https://academic.oup.com/bioinformatics/article/28/19/2520/290322) on its history.

Both Snakemake and Nextflow are using dataflow programming language. Think of this as mapping on inputs and outputs; jobs that have gathered all required inputs will immediately run. Print out the directed acyclic graph will help you visualize your jobs.

Config
```console
$ cat config.yaml
samples:
  - Sample1
  - Sample2
  - Sample3

reference: /data/home/misc/Snakemake_vs_nextflow/code_club/snakemake_version2/input_data/reference/ARS-UCD2.0_demo.fa

$ cat environment.yaml
channels:
  - bioconda
  - conda-forge
dependencies:
  - snakemake-minimal >=7.3
  - jinja2
  - matplotlib
  - graphviz
  - bcftools =1.15
  - samtools =1.15
  - bwa =0.7.17
  - pysam =0.19
  - pygments
```

How to run snakemake?
```console
conda activate lloydcondaenv
snakemake -n -s Snakefile --configfile config.yaml
snakemake -n -s Snakefile --configfile config.yaml --dag | dot -Tsvg > snakedag.svg
snakemake --cores 16 -s Snakefile --configfile config.yaml
```


The snake dag ... ![snake dag](snakedag.svg)

The snakefile
```python3
rule all:
	input:
		"results/multiqc/multiqc_report.html",
		#expand("results/trim_galore/{sample}_{read}_val_{read}.fq.gz", sample=config["samples"], read=["1", "2"]),
		#expand("results/bwa/{sample}.bam", sample=config["samples"]),
		expand("results/bcftools/{sample}.vcf.gz", sample=config["samples"])

rule raw_fastqc:
	input:
		r1="/data/home/misc/Snakemake_vs_nextflow/code_club/snakemake_version2/input_data/raw_fastq/{sample}.R1.fastq.gz",
		r2="/data/home/misc/Snakemake_vs_nextflow/code_club/snakemake_version2/input_data/raw_fastq/{sample}.R2.fastq.gz"

	output:
		"results/fastqc/raw/{sample}.R1_fastqc.html",
		"results/fastqc/raw/{sample}.R2_fastqc.html"

	params:
		threads=4

	log:
		"logs/fastqc/raw/{sample}.log"

	shell:
		"""
		fastqc -t {params.threads} {input.r1} {input.r2} --outdir results/fastqc/raw &> {log}
		"""

rule trim_galore:
	input:
		r1="/data/home/misc/Snakemake_vs_nextflow/code_club/snakemake_version2/input_data/raw_fastq/{sample}.R1.fastq.gz",
		r2="/data/home/misc/Snakemake_vs_nextflow/code_club/snakemake_version2/input_data/raw_fastq/{sample}.R2.fastq.gz"

	output:
		"results/trim_galore/{sample}.R1_val_1.fq.gz",
		"results/trim_galore/{sample}.R2_val_2.fq.gz"

	params:
		threads=4

	log:
		"logs/trim_galore/{sample}.log"

	shell:
		"""
		trim_galore --quality 30 \
    	--length 75 \
    	--clip_R1 1 \
    	--clip_R2 1 \
    	--three_prime_clip_R1 1 \
    	--three_prime_clip_R2 1 \
    	--paired \
    	--cores {params.threads} \
		{input.r1} {input.r2} -o results/trim_galore &> {log}
		"""

rule trimmed_fastqc:
	input:
		r1="results/trim_galore/{sample}.R1_val_1.fq.gz",
		r2="results/trim_galore/{sample}.R2_val_2.fq.gz"

	output:
		"results/fastqc/trimmed/{sample}.R1_val_1_fastqc.html",
		"results/fastqc/trimmed/{sample}.R2_val_2_fastqc.html"

	params:
		threads=4

	log:
		"logs/fastqc/trimmed/{sample}.log"

	shell:
		"""
		fastqc -t {params.threads} {input.r1} {input.r2} --outdir results/fastqc/trimmed &> {log}
		"""

rule bwa:
	input:
		r1="results/trim_galore/{sample}.R1_val_1.fq.gz",
		r2="results/trim_galore/{sample}.R2_val_2.fq.gz"

	output:
		"results/bwa/{sample}.bam"

	params:
		threads=4

	log:
		"logs/bwa/{sample}.log"

	shell:
		"""
		bwa mem -t {params.threads} \
		-R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:{wildcards.sample}" \
		{config[reference]} \
		{input.r1} {input.r2} | \
		samtools view -b - | samtools sort -o {output} &> {log}
		samtools index {output}
		"""

rule markDuplicates:
	input:
		"results/bwa/{sample}.bam"

	output:
		"results/picard/{sample}.markDup.bam"

	log:
		"logs/picard/{sample}.log"

	shell:
		"""
		java -jar -Xmx3g /home/lloyd/Software/picard/build/libs/picard.jar MarkDuplicates  \
        I={input} \
        O=results/picard/{wildcards.sample}.markDup.bam \
        M=results/picard/{wildcards.sample}.markDup.metrics.txt \
        VALIDATION_STRINGENCY=LENIENT \
        REMOVE_DUPLICATES=false \
        OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 \
        CREATE_INDEX=false
    samtools index results/picard/{wildcards.sample}.markDup.bam
	"""

rule callVariants:
	input:
		"results/picard/{sample}.markDup.bam"

	output:
		"results/bcftools/{sample}.vcf.gz"

	threads: 4

	log:
		"logs/bcftools/{sample}.log"

	shell:
		"""
		bcftools mpileup -f {config[reference]} {input} | \
    	bcftools call --threads {threads} -m --output-type z | \
    	bcftools filter --threads {threads} --output-type z -s LowQual -e 'QUAL<20 || DP<100' > results/bcftools/{wildcards.sample}.vcf.gz
    	bcftools index results/bcftools/{wildcards.sample}.vcf.gz
		bcftools stats results/bcftools/{wildcards.sample}.vcf.gz > results/bcftools/{wildcards.sample}.vcf.stats
		"""

rule multiqc:
	input:
		expand("results/fastqc/raw/{sample}.R{read}_fastqc.html", sample=config["samples"], read=["1", "2"]),
		expand("results/fastqc/trimmed/{sample}.R{read}_val_{read}_fastqc.html", sample=config["samples"], read=["1", "2"]),
		#expand("results/trim_galore/{sample}_{read}_val_{read}.fq.gz", sample=config["samples"], read=["1", "2"]),

	output:
		"results/multiqc/multiqc_report.html"

	shell:
		"""
		multiqc results/fastqc/raw results/fastqc/trimmed -o results/multiqc
		"""
```

After run directory structure

```console
$ ls -lah
total 70K
drwxrwsr-x  6 lloyd lloyd   11 Dec 14 11:39 .
drwxrwsr-x 10 lloyd lloyd   10 Dec 13 18:52 ..
-rw-rw-r--  1 lloyd lloyd  164 Dec 13 22:46 config.yaml
-rw-rw-r--  1 lloyd lloyd  203 Dec 13 18:52 environment.yaml
drwxrwsr-x  4 lloyd lloyd    4 Dec 13 22:41 input_data
drwxrwsr-x  7 lloyd lloyd    7 Dec 14 11:43 logs
-rw-rw-r--  1 lloyd lloyd  218 Dec 13 23:08 mylog
drwxrwsr-x  8 lloyd lloyd    8 Dec 14 11:43 results
-rw-rw-r--  1 lloyd lloyd  17K Dec 14 11:31 snakedag.svg
-rw-rw-r--  1 lloyd lloyd 4.1K Dec 13 23:17 Snakefile
drwxrwsr-x 11 lloyd lloyd   11 Dec 14 11:39 .snakemake
```

<a name="nextflow"></a>
## Nextflow

![Nextflow](nextflow.png)

I recommend [Nextflow tutorial](https://www.nextflow.io/docs/latest/getstarted.html) and this simple example with [blast](https://www.nextflow.io/example3.html).

```console
$ cat insyb2023.yaml
name: insyb2023
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - bcftools
  - bwa
  - cutadapt
  - fastqc=0.12.1=hdfd78af_0
  - sambamba
  - samtools=1.10
  - multiqc
  - trim-galore

mamba env create -f insyb2023.yaml
```

<a name="Acknowledgements"></a>
## Acknowledgements

I thank the Davies coding club members for practicing on code similar to the one presented here. In particular, special thanks to Callum MacPhillamy, my post-doc, who coded the majority of the pipelines presented here.
