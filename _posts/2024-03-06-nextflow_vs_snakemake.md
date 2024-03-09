--
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

During the last 2 years, I have come across Nextflow and Snakemake, which both seem to have the aim to provide an infrastructure for "reproducible and scalable" workflows. In early 2023 I started using Snakemake, and in practice it effectively seemed to provide a cleaner and more organized way to develop and run pipelines. In particular, if you need to run different tools that need their own Conda environment, with Snakemake you can seamlessly define such environments and it will automatically switch them accordingly at execution time. This, however, also seems to be one of its drawbacks as it needs Conda to do so, and this can take ages some times (the always there "Solving environment".. Just type "conda so" in Google and you will see that the most searched phrase is related to this). Though there is a way to make it run with Mamba, which is faster, it can get confusing if you haven't work with Python and Conda environments before. On the other hand, Nextflow, as we will see later, seems simpler, and it doesn't require jumping into the whole snakes universe of Conda-Mamba-Snakemake.

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



