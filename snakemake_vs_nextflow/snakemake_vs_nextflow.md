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

### Reference genome

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
