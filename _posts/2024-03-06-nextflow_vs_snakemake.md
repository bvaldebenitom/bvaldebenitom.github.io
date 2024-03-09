---
layout: post
title:  Nextflow vs Snakemake
date: 2024-03-06 05:55:00
description: A comparison of Nextflow and Snakemake in the context of RNA-Seq analysis
tags: workflows rnaseq
categories: reviews
thumbnail: assets/img/02_nextflow_vs_snakemake/figure01.png
toc:
  beginning: true
---

## Introduction
Throughout my career, I have always enjoyed the development of pipelines as a mean to efficiently process hundreds of files efficiently. I usually do this directly on Bash, as these pipelines connect the inputs and outputs of different software. I have been able to publish some of these tools in several journals, and at this point I'm always wondering whether they will work in all environments for all users. Although I have made a conscious effort to make them reproducible, they still fail at some point. Another issue is with Bash pipelines, is that I always have to make blocks of "if-else" statements just to check whether the inputs and outputs were created, which is very cumbersome.

During the last 2 years, I have come across Nextflow and Snakemake, which both seem to have the aim to provide an infrastructure for "reproducible and scalable" workflows. In early 2023 I started using Snakemake, and in practice it effectively seemed to provide a cleaner and more organized way to develop and run pipelines. In particular, if you need to run different tools that need their own Conda environment, with Snakemake you can seamlessly define such environments and it will automatically switch them accordingly at execution time. This, however, also seems to be one of its drawbacks as it needs Conda to do so, and this can take ages some times (the always there "Solving environment".. Just type "conda so" in Google and you will see that the most searched phrase is related to this). Though there is a way to make it run with Mamba, which is faster, it can get confusing if you haven't worked with Python and Conda environments before. On the other hand, Nextflow, as we will see later, seems simpler, and it doesn't require jumping into the whole snakes universe of Conda-Mamba-Snakemake.

To get more familiarized with both Nextflow and Snakemake, I adapted a simplified version of my RNA-Seq workflow from scratch to each of these frameworks. In this post I will be discussing some of my findings and what I would recommend for someone wanting to implement them. 

## The workflow

For this workflow, I simulated reads in FASTQ formats: 3 paired-end libraries and 1 single-end library. Then, these reads would need to be quality checked with FastQC, aligned against the human genome with STAR, and finally, processed with Telescope to get expression estimates for Transposable Elements (Figure 1).

<div class="row mt-3">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.html path="assets/img/02_nextflow_vs_snakemake/figure01.png" class="img-fluid rounded z-depth-1" zoomable=true %}
    </div>
</div>

There are 2 considerations for the pipeline:
1. It should automatically detect whether the libraries are single-end or paired-end.
2. Telescope needs to be run on its own Conda environment (I have already built it), so it should be able to reuse the same existing environment.

I will show the implementation first on Snakemake and latter on Nextflow.

### Snakemake implementation

#### Setup

The basic commands for [installing Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) are:

````markdown
conda install -n base -c conda-forge mamba
mamba create -c conda-forge -c bioconda -n snakemake snakemake
mamba activate snakemake
snakemake --help
````

So, right off the bat, we already need to have Conda installed, and before installing Snakemake, they recommend to install Mamba **using Conda**. If you do this, you might lose a good amount of time. I recommend installing Mamba or Micromamba following their instructions. In my case, I chose to do this with Micromamba, which I started using thanks to a good friend. This is how the commands look:

````markdown
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
micromamba create -c conda-forge -c bioconda -n snakemake snakemake
micromamba activate snakemake
````

#### Writing the workflow

First, we will define paths to the FastQC and STAR binaries and other parameters:

````markdown
params_fastqc = "Snakemake_Nextflow/FastQC/fastqc"
params_fastqc_memory = "10000"
params_fastqc_threads = 1

params_star_index = "STAR_GenomeIndex/hg38_STAR_2.7.11b"
params_star = "STAR_2.7.11b/Linux_x86_64/STAR"
params_star_threads = 21

params_telescope_gtf = "hg38_all.gtf"
````

Next, we will set some variables to guide the different steps of the pipeline. Snakemake allows the use of Python commands, so I'm taking advantage of this. First, the "glob" library is imported, and the "samples" variable is created by listing all files with extension "fq". Concurrently, through the use of a regular expression, we substitute the trailing _1 or _2 of the paired-end files. In the regular expression we add the quantifier "{0,1}" which makes it also match single-end files that just end on ".fq", without a trailing _1 or _2. 

````markdown
import glob

samples = list(set([re.sub("(_[12]){0,1}.fq","",x) for x in glob.glob("*.fq")]))
# ['hg38_tes_random_n1000_simulated_pe_reads_l150_f15_SE', 'hg38_tes_random_n1000_simulated_pe_reads_l150_f20_m1000_s10', 'hg38_tes_random_n1000_simulated_pe_reads_l150_f20_m1100_s10', 'hg38_tes_random_n1000_simulated_pe_reads_l150_f20_m900_s10']
````

Similarly, the "fastq_basenames" variable is created:
````markdown
fastq_basenames = [re.sub(".fq","",x) for x in glob.glob("*.fq")]
# ['hg38_tes_random_n1000_simulated_pe_reads_l150_f20_m1100_s10_1', 'hg38_tes_random_n1000_simulated_pe_reads_l150_f20_m1000_s10_2', 'hg38_tes_random_n1000_simulated_pe_reads_l150_f20_m1100_s10_2', 'hg38_tes_random_n1000_simulated_pe_reads_l150_f15_SE', 'hg38_tes_random_n1000_simulated_pe_reads_l150_f20_m900_s10_2', 'hg38_tes_random_n1000_simulated_pe_reads_l150_f20_m1000_s10_1', 'hg38_tes_random_n1000_simulated_pe_reads_l150_f20_m900_s10_1']
````

The variables "pe_samples" and "se_samples" follow a similar logic, with the exception that the first one will only contain the samples that are paired-end and the second one the ones that are single-end. These variables are used in the function "input_fastqs" which will process a Snakemake input and identify which type is the sample.

````markdown
pe_samples = list(set([re.sub("(_[12]){0,1}.fq","",x) for x in glob.glob("*_[12].fq")]))
# ['hg38_tes_random_n1000_simulated_pe_reads_l150_f20_m1000_s10', 'hg38_tes_random_n1000_simulated_pe_reads_l150_f20_m1100_s10', 'hg38_tes_random_n1000_simulated_pe_reads_l150_f20_m900_s10']

se_samples = [re.sub(".fq","",x) for x in glob.glob("*.fq") if re.search("[^_12].fq",x)]
# ['hg38_tes_random_n1000_simulated_pe_reads_l150_f15_SE']

def input_fastqs(wildcards):
        print(wildcards.sample)
        if wildcards.sample in se_samples:
                return(f'{wildcards.sample}.fq')
        if wildcards.sample in pe_samples:
                return([f'{wildcards.sample}_1.fq',f'{wildcards.sample}_2.fq'])
````

















