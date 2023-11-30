# Welcome to THE_Aligner

THE_Aligner is a Snakemake pipeline developed by Isaac Vock. It is designed to process fastq files from various flavors of RNA-seq experiments and align them to a reference genome.

## What THE_Aligner does

The pipeline includes the following steps:

1. Trim adapters with [fastp](https://github.com/OpenGene/fastp)
    * Fastqs will also be unzipped with [pigz]() if gzipped. If this is the case, the unzipped fastqs are temporary files that get removed once the pipeline steps using them have finished running. This saves on disk space
1. Assess fastqs with [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
1. Align fastqs
    * Includes the strictly splice unaware aligners [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) or[bowtie2](https://github.com/BenLangmead/bowtie2).
    * Also includes the splice aware aligners [star](https://github.com/alexdobin/STAR) and [hisat2](https://github.com/DaehwanKimLab/hisat2).
    * Alignment statistics are generated with [bamtools](https://github.com/pezmaster31/bamtools).
    * Alignment indices can also be automatically built for all implemented aligners.
1. Sort bam files with [samtools](http://www.htslib.org/doc/samtools-sort.html)
1. Generate coverage files
    * Bedgraph files created with [bedtools](https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html)
    * BigWig files created with [bedGraphtoBigWig](https://www.encodeproject.org/software/bedgraphtobigwig/)

In addition, if quantification of annotated transcripts is all you are interested in, THE_Aligner also implements the two most popular pseudo aligners, [kallisto](https://pachterlab.github.io/kallisto/about) and [salmon](https://combine-lab.github.io/salmon/). In this case, bam files are not generated so the sorting and coverage file creation steps are skipped. Once again, indices for these pseudoaligners can be automatically generated.


## Requirements for THE_Aligner

THE_Aligner uses the workflow manager [Snakemake](https://snakemake.readthedocs.io/en/stable/). The minimal version of Snakemake is techncially compatible with Windows, Mac, and Linux OS, but several of the software dependencies are only Mac and Linux compatible. If you are a Windows user like me, don't sweat it, I would suggest looking to the Windows subsystem for linux which can be easily installed (assuming you are running Windows 10 version 2004 or higher).

In addition, you will need Git installed on your system so that you can clone this repository. Head to [this link](https://git-scm.com/downloads) for installation instructions if you don't already have Git.

## Getting Started

There are several ways to run THE_Aligner:

1. [Deploying with Snakedeploy](deploy.md) (recommended)
    * A Simon Lab/Yale specific version of these instructions are [here](simon.md). While highly specific, this also includes instructions for optimized deployment on a cluster using a Slurm scheduler, and could thus be of some general use.
1. [Cloning the repo locally](alt.md)