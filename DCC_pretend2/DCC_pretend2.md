# Davies coding club
---
*20231109*

## Table of Contents
* [DCC_pretend](#dcc)

<a name="dcc"></a>
## DCC_pretend
Pretend making snakemake pipeline together

Useful snakemake related commands

```console
screen
conda activate snakemake-tutorial
snakemake -n -s Snakefile --configfile config.yaml
snakemake -n -s Snakefile --configfile config.yaml --forceall
snakemake --cores 4 -s Snakefile --configfile config.yaml
snakemake --cores 4 -s Snakefile --configfile config.yaml --forceall
snakemake --cores 16 -s Snakefile --configfile config.yaml --forceall
snakemake --cores 4 -p markdup_results/Sample?/Sample?_marked.bam -s Snakefile --configfile config.yaml --force
```

# Rules for step 1 (mapping fastq to ref)

```console
rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/A.fastq"
    output:
        "mapped_reads/A.bam"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"
```

# Rules for step 3 (Indexing read alignments and visualizing the DAG)

```console
rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"
        
 ```

 # Rule for step 5 (Calling genomic variants) 
 ```console

 SAMPLES = ["A", "B"]
 rule bcftools_call:
    input:
        fa="data/genome.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=SAMPLES),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES)
    output:
        "calls/all.vcf"
    shell:
        "bcftools mpileup -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"
 
 ```      
        
        
        
        
        
        
        
        