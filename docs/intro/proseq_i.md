# Welcome to PROseq_etal

PROseq_etal is a Snakemake pipeline developed by Isaac Vock. It is designed to process fastq files from PRO-seq and ChIP-seq experiments. Future development will see it also accomodate a number of other related modalities (e.g., NET-seq, ATAC-seq, etc.).

## What PROseq_etal does

The pipeline includes the following steps:

1. Trim adapters with [fastp](https://github.com/OpenGene/fastp)
    * Fastqs will also be unzipped with [pigz]() if gzipped. If this is the case, the unzipped fastqs are temporary files that get removed once the pipeline steps using them have finished running. This saves on disk space
1. Assess fastqs with [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
1. Align fastqs
    * Can use [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) or [bowtie2](https://github.com/BenLangmead/bowtie2)
    * Alignment indices can also be built automatically if not provided by user
    * Alignment statistics are generated with [bamtools](https://github.com/pezmaster31/bamtools)
1. Sort bam files with [samtools](http://www.htslib.org/doc/samtools-sort.html)
1. Generate coverage files
    * Bedgraph files created with [bedtools](https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html)
    * BigWig files created with [bedGraphtoBigWig](https://www.encodeproject.org/software/bedgraphtobigwig/)
1. Call peaks with [MACS2](https://github.com/macs3-project/MACS/tree/master)
    * If processing ChIP-seq data, fold-enrichment relative to Input is also computed with MACS2
    * In addition, for ChIP-seq data, fold-enrichment bedGraphs are converted to bigWigs
1. Identify and annotate peaks/transcripts with [HOMER](http://homer.ucsd.edu/homer/)
1. Quantify gene body and pause site coverage with [HTSeq](https://htseq.readthedocs.io/en/master/htseqcount.html) (PRO-seq DATA ONLY)
    * Pause indices are also calculated with a custom R script

## Requirements for PROseq_etal

PROseq_etal uses the workflow manager [Snakemake](https://snakemake.readthedocs.io/en/stable/). The minimal version of Snakemake is techncially compatible with Windows, Mac, and Linux OS, but several of the software dependencies are only Mac and Linux compatible. If you are a Windows user like me, don't sweat it, I would suggest looking to the Windows subsystem for linux which can be easily installed (assuming you are running Windows 10 version 2004 or higher).

In addition, you will need Git installed on your system so that you can clone this repository. Head to [this link](https://git-scm.com/downloads) for installation instructions if you don't already have Git.

## Getting Started

There are several ways to run PROseq_etal:

1. [Deploying with Snakedeploy](deploy.md) (recommended)
    * A Simon Lab/Yale specific version of these instructions are [here](simon.md). While highly specific, this also includes instructions for optimized deployment on a cluster using a Slurm scheduler, and could thus be of some general use.
    * Information about configuring PROseq_etal is [here](../configuration/proseq_c.md).
1. [Cloning the repo locally](alt.md)