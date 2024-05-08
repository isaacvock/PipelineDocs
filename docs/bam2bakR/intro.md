# Welcome to bam2bakR

[bam2bakR](https://github.com/simonlabcode/bam2bakR) is a Snakemake implementation of a portion of the [TimeLapse pipeline](https://bitbucket.org/mattsimon9/timelapse_pipeline/src/master/) developed by the [Simon lab](https://simonlab.yale.edu/) at Yale. The contributors to the original pipeline are Matthew Simon, Jeremy Schofield, Martin Machyna, Lea Kiefer, and Joshua Zimmer.

## What bam2bakR does

bam2bakR includes the basic TimeLapse pipeline functionality downstream of alignment. Thus, the input to bam2bakR is a set of .bam files and the output is a cB.csv file, which the following columns:

* sample - Sample name
* rname - Chromosome name
* sj - Logical: TRUE if read overlaps an exon-exon splice junction
* XF - Exonic feature: Gene ID if the read mapped solely to exonic parts of a gene
* GF - Gene feature: Gene ID when read is aligned to any part of a gene (intronic or exonic)
* nT - Number of uridines in RNA from which read was derived (additional columns are added if considering other types of mutations, like G-to-A mutations)
* TC - Number of called T-to-C mutations in read (additional columns are added for other types of mutations)
* n - Number of reads which have the identical set of values described above

Additional columns can be included and any of these columns can be omitted per your choice.

## Requirements for bam2bakR

bam2bakR uses the workflow manager [Snakemake](https://snakemake.readthedocs.io/en/stable/). The minimal version of Snakemake is techncially compatible with Windows, Mac, and Linux OS, but several of the software dependencies are only Mac and Linux compatible. If you are a Windows user like me, don't sweat it, I would suggest looking to the Windows subsystem for linux which can be easily installed (assuming you are running Windows 10 version 2004 or higher). 

To install Snakemake, I would suggest installing the package manager [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) (e.g., by installing [miniforge](https://github.com/conda-forge/miniforge)) and using that to create an environment containing Snakemake and [Snakedeploy](https://github.com/snakemake/snakedeploy). Snakedeploy is the tool needed for the preferred method of deploying bam2bakR on your system. Once mamba is installed on your system, this can be accomplished in one line: `mamba create -c conda-forge -c bioconda --name deploy_snakemake snakemake snakedeploy`, where I have called the new environment "deploy_snakemake". Run `conda activate deploy_snakemake` and you can now use Snakemake and Snakedeploy!

In addition, you will need Git installed on your system so that you/Snakedeploy can clone this repository. Head to [this link](https://git-scm.com/downloads) for installation instructions if you don't already have Git.

Finally, bam2bakR requires bam files from stranded library preps (i.e., information about whether a given read represents the original sequence of the RNA or its reverse complement must be retained). In addition, make sure that the aligner you used was configured to output bam files with an MD tag, which will keep track of the reference nucleotide at any mismatches between the aligned read and the reference genome. While this tag is not always included by default, most aligners have an option to add this tag to the list of default tags.


## Getting Started

There are a number of ways to use bam2bakR:

1. [Deploying with Snakedeploy](../deploy.md) (recommended)
    * A Simon Lab/Yale specific version of these instructions are [here](../simon.md). While highly specific, this also includes instructions for optimized deployment on a cluster using a Slurm scheduler, and could thus be of some general use.
    * Information about configuring bam2bakR is [here](../bam2bakR/configuration.md).
1. [Cloning the repo locally](../alt.md)