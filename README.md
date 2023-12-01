# PipelineDocs
One-stop-shop for documentation of my Snakemake pipelines.

Link: [https://pipelinedocs.readthedocs.io/en/latest/](https://pipelinedocs.readthedocs.io/en/latest/)

Pipelines currently discussed in the documentation:

1. [bam2bakR](https://github.com/simonlabcode/bam2bakR/blob/main/README.md)
    - Snakemake reimplementation of a pipeline the Simon lab developed to process TimeLapse-seq and similar data.
2. [THE_Aligner](https://github.com/isaacvock/THE_Aligner)
    - Pipeline for aligning almost any kind of sequencing data.
    - Implements several aligners (bowtie2, bwa-mem2, STAR, HISAT2).
    - Also implements popular pseudoaligners for rapid quantification of transcript/gene abundances without alignment (kallisto and salmon).
3. [PROseq_etal](https://github.com/isaacvock/PROseq_etal)
    - Pipeline for processing PRO-seq and ChIP-seq data. Originally based off of ChIP-seq pipeline developed by the Simon lab.
