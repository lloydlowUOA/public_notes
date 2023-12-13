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

$ ls -lah raw_fastq/
total 597M
drwxrwsr-x  2 lloyd lloyd    9 Dec  5 14:52 .
drwxrwsr-x 10 lloyd lloyd   10 Dec 13 18:52 ..
-rw-rw-r--  1 lloyd lloyd 6.1K Dec  5 14:29 .DS_Store
-rw-rw-r--  1 lloyd lloyd  82M Dec  5 14:29 Sample1.R1.fastq.gz
-rw-rw-r--  1 lloyd lloyd  83M Dec  5 14:29 Sample1.R2.fastq.gz
-rw-rw-r--  1 lloyd lloyd 127M Dec  5 14:29 Sample2.R1.fastq.gz
-rw-rw-r--  1 lloyd lloyd 128M Dec  5 14:29 Sample2.R2.fastq.gz
-rw-rw-r--  1 lloyd lloyd  89M Dec  5 14:29 Sample3.R1.fastq.gz
-rw-rw-r--  1 lloyd lloyd  90M Dec  5 14:29 Sample3.R2.fastq.gz

$ ls -lah reference/
total 192M
drwxrwsr-x  2 lloyd lloyd   9 Dec  5 14:29 .
drwxrwsr-x 10 lloyd lloyd  10 Dec 13 18:52 ..
-rw-rw-r--  1 lloyd lloyd 95M Dec  5 14:29 ARS-UCD2.0_demo.fa
-rw-rw-r--  1 lloyd lloyd  97 Dec  5 14:29 ARS-UCD2.0_demo.fa.amb
-rw-rw-r--  1 lloyd lloyd  77 Dec  5 14:29 ARS-UCD2.0_demo.fa.ann
-rw-rw-r--  1 lloyd lloyd 93M Dec  5 14:28 ARS-UCD2.0_demo.fa.bwt
-rw-rw-r--  1 lloyd lloyd  53 Dec  5 14:28 ARS-UCD2.0_demo.fa.fai
-rw-rw-r--  1 lloyd lloyd 24M Dec  5 14:29 ARS-UCD2.0_demo.fa.pac
-rw-rw-r--  1 lloyd lloyd 47M Dec  5 14:28 ARS-UCD2.0_demo.fa.sa

</p>
</details>

<a name="shell"></a>
## Shell scripting



<a name="snakemake"></a>
## Snakemake

![Snakemake](snakemake.png)

I recommend [Snakemake tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html)

<a name="nextflow"></a>
## Nextflow

![Nextflow](nextflow.png)

I recommend [Nextflow tutorial](https://www.nextflow.io/docs/latest/getstarted.html)

<a name="Acknowledgements"></a>
## Acknowledgements

I thank the Davies coding club members for practicing on code similar to the one presented here. In particular, special thanks to Callum MacPhillamy, my post-doc, who coded the majority of the pipelines presented here.
