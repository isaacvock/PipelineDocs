## Configuring bam2bakR

In the `config/` directory you will find a file named `config.yaml`. If you open it in a text editor, you will see several parameters which you can alter to your heart's content. The first parameter that you have to set is at the top of the file:

``` yaml
samples:
  WT_1: data/samples/WT_replicate_1.bam
  WT_2: data/samples/WT_replicate_2.bam
  WT_ctl: data/samples/WT_nos4U.bam
  KO_1: data/samples/KO_replicate_1.bam
  KO_2: data/samples/KO_replicate_2.bam
  KO_ctl: data/samples/KO_nos4U.bam
```
`samples` is the list of sample IDs and paths to .bam files that you want to process. Delete the existing sample names and paths and add yours. The sample names in this example are `WT_1`, `WT_2`, `WT_ctl`, `KO_1`, `KO_2`, and `KO_ctl`. These are the sample names that will show up in the `sample` column of the output cB.csv file. The `:` is necessary to distinguish the sample name from what follows, the path to the relevant bam file. Note, the path is NOT an absolute path, it is relative to the directory that you deployed to (i.e., `workdir` in this example). Thus, in this example, the bam files are located in a directory called `samples` that is inside of a directory called `data` located in `workdir`. Your data can be wherever you want it to be, but it might be easiest if you put it in a `data` directory inside the bam2bakR directory as in this example. 

As another example, imagine that the `data` directory was in the directory that contains `workdir`, and that there was no `samples` subdirectory inside of `data`. In that case, the paths would look something like this:

``` yaml
samples:
  WT_1: ../data/WT_replicate_1.bam
  WT_2: ../data/WT_replicate_2.bam
  WT_ctl: ../data/WT_nos4U.bam
  KO_1: ../data/KO_replicate_1.bam
  KO_2: ../data/KO_replicate_2.bam
  KO_ctl: ../data/KO_nos4U.bam
```
where `../` means navigate up one directory. 

The next parameter you have to set denotes the sample names of any -s4U control samples (i.e., samples that were not fed s4U or a similar metabolic label):

``` yaml
control_samples: ["WT_ctl", "KO_ctl"]
```

In this case, the samples named WT_ctl and KO_ctl are the -s4U control samples. -s4U controls will be used to call any single nucleotide polymorphisms (SNPs) in your cell line so as to avoid conflating them with T-to-C mutations induced by the nucleotide recoding chemistry. 

The third crucial parmaeter immediately follows:

``` yaml
annotation: data/annotation/GRCh38.gtf
```
This is the path to the GTF annotation file for the genome that reads were mapped to. The same rules apply when it comes to specifying this path.

Finally, the path to the genome fasta file that you used for alignment must also be specified:

``` yaml
genome_fasta: data/genome/GRCh38.fasta
```

The other parameters that can be altered are:

* `strandedness`: whether the first read in a pair (or the only read if single-end) represents the original sequence of the RNA (F), or its reverse complement (R). For example, set this parameter to "F" if your library is an FR paired-end library, and "R" if it is an RF paired-end library.
* `FORMAT`: whether the reads are paired-end (PE) or single-end (SE).
* `mut_tracks`: the type of mutation (e.g., T-to-C mutations) that sequencing tracks will be colored by. If you are most interested in the T-to-C mutational content, then `mut_tracks` should be TC. If G-to-A, then `mut_tracks` should be GA. If both, then `mut_tracks` should be "TC,GA".
* `minqual`: Minimum base quality to call it a mutation.
* `spikename`: If spike-ins are present, this should be a string that is common to all gene_ids for spike-in transcripts in annotation gtf. For example, in Ensembl annotations for Drosophila melanogaster, all gene_ids start with "FBgn". Therefore, if you have Drosophila spike-ins, `spikename` should be "FBgn".
* `normalize`: If True, then scale factor calculated with edgeR is used to normalize sequencing tracks.
* `WSL`: whether you are running this on the Windows subsystem for linux (0 = yes; 1= no)

 
Edit the values in the config file as necessary and move on to the last step.


## bam2bakR with fastq input

Despite its name, bam2bakR is also able to take fastq files as input, trimming and aligning them before running them through the standard bam2bakR pipeline. The current implementation of this is a bit rough around the edges and inflexible though. For aligning fastq files to create bam files that can be passed to bam2bakR, I will shamelessly suggest my own Snakemake pipeline called [THE_Aligner](https://github.com/isaacvock/THE_Aligner) which implements a number of different aligners, makes use of automated unit testing, is more configurable, and includes QC of fastqs and output bam files that bam2bakR does not currently implement. A new and improved version of the full fastq-to-cB file pipeline implemented in bam2bakR is currently under development at [fastq2EZbakR](https://github.com/isaacvock/fastq2EZbakR), and is almost fully operational (as of 2/1/2024; pipeline is fully implemented, though slowly migrating some of its functionality to bam2bakR). Documentation for configuring bam2bakR's fastq processing strategy is included here for posterity's sake.

All of the parameters discussed for the case of bam file input apply equally well to fastq inputs; the only difference is that the paths in the samples parameter are interpreted a bit differently:

``` yaml
samples:
  WT_1: data/fastq/WT_1
  WT_2: data/fastq/WT_2
  WT_ctl: data/fastq/WT_ctl
  KO_1: data/fastq/KO_1
  KO_2: data/fastq/KO_2
  KO_ctl: data/fastq/KO_ctl
```
In this case, `samples` is still the list of sample IDs, but the paths are to directories containing .fastq files that you want to process, rather than to individual files. **NOTE: each directory must contain a single fastq file or pair of fastq files**.   

Configurable parameters that uniquely impact the fastq processing are:

* `HISAT2`: Path to directory containing hisat2 indices. These must be pre-built by the user. See [here](https://daehwankimlab.github.io/hisat2/manual/) for details.
* `HISAT_3N`: Path to directory containing hisat-3n indices. These must be pre-built by the user. See [here](https://daehwankimlab.github.io/hisat2/hisat-3n/) for details.
* `STAR_index`: Path to directory containing STAR indices. These can be created by the pipeline if `build_star` is set to True. They will be created at the directory specified by `STAR_index`.
* `use_hisat3n`: If True, HISAT-3N will be used for alignment.
* `use_star`: If True, and `use_hisat3n` is False, STAR will be used for alignment. If both `use_hisat3n` and `use_star` are false, then HISAT2 will be used.
* `build_star`: If True, STAR indices will be built automatically at path specified in `STAR_index` if not provided by user.
* `hisat3n_path`: Path to HISAT-3N executable; HISAT-3N cannot be installed automatically like all other dependencies as it is not installable via conda.
* `chr_tag`: If True, chr is added to chromosome names during alignment (HISAT2 and HISAT-3N only). Useful when aligner index is number-based but annotation file has "chr" in its seqnames.
* `Yale`: If True, HISAT-3N can be loaded as a module available on the McCleary cluster. Otherwise, HISAT-3N must be installed manaully by the user and the path to the exectuable provided.
* `flattened`: This parameter should be kept to False unless you know what you are doing.
* `adapter`: Code snippet that is passed to [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) for adapter trimming.
* `cutadapt_extra`: Extra, optional parameters passed to Cutadapt.
* `star_extra`: Extra, optional parameters passed to STAR for alignment.
* `hisat2_extra`: Extra, optional parameters passed to HISAT2 for alignment.