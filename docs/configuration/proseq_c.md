## Configuring PROseq_etal

In the `config/` directory you will find a file named `config.yaml`. If you open it in a text editor, you will see several parameters which you can alter to your heart's content. The first parameter that you have to set is at the top of the file:

``` yaml
samples:
  WT_1: data/fastq/WT_1
  WT_2: data/fastq/WT_2
  WT_ctl: data/fastq/WT_ctl
  KO_1: data/fastq/KO_1
  KO_2: data/fastq/KO_2
  KO_ctl: data/fastq/KO_ctl
```
`samples` is the list of sample IDs and paths to .bam files that you want to process. Delete the existing sample names and paths and add yours. The sample names in this example are `WT_1`, `WT_2`, `WT_ctl`, `KO_1`, `KO_2`, and `KO_ctl`. These are the sample names that will append many of the files output by PROseq_etal. The `:` is necessary to distinguish the sample name from what follows, the path to the relevant bam file. Note, the path can be absolute (e.g., ~/path/to/fastqs/) or relative to the directory that you deployed to (i.e., `workdir` in this example). In the example config, the paths specified are relative. Thus, in this example, the bam files are located in a directory called `samples` that is inside of a directory called `data` located in `workdir`. Your data can be wherever you want it to be, but it might be easiest if you put it in a `data` directory inside the PROseq_etal directory as in this example. 

As another example, imagine that the `data` directory was in the directory that contains `workdir`, and that there was no `samples` subdirectory inside of `data`. In that case, the relative paths would look something like this:

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

The next parameter you have to set denotes the experimental method used:

``` yaml
method: "ChIPseq"
```

The current options are "ChIPseq" and "PROseq", with more to come! 

The third parameter is only relevant if you are analyzing ChIPseq data, and it identifies the relevant Input control samples for each enrichment sample:

``` yaml
controls:
  WT_1: WT_ctl
  WT_2: WT_ctl
  KO_1: KO_ctl
  KO_2: KO_ctl
```
The "keys" (what is to the left of the ":") are sample IDs from `samples:`. The sample IDs in this section should only correspond to the IDs for enriched samples. The "values" (what is to the right of the ":") is the sample ID for the relevant Input sample. In this example, the sample labeled WT_ctl is the Input sample from which the WT_1 enrichment sample was derived. Fold enrichment tracks for WT_1 will be calculated using WT_ctl as the Input reference.

The remaining somewhat more self-explanatory required parameters are:

* `genome`: Path to genome fasta file to be used for alignment.
* `aligner`: Determines which aligner will be used (options are "bwa-mem2" and "bowtie2").
* `indices`: Path to aligner indices. These will be created at this path if not already present.
* `annotation`: Path to annotation gtf file to be used by HTSeq.
* `PI_gtf`: Path to pause site annotation gtf file to be used to calculate pause indices with HTSeq and custom scripts. This will be created automatically at this path if it does not already exist.
* `strandedness`: HTSeq parameter specifying library strandedness. Options are "reverse", "yes", or "no". See config comments and/or [HTSeq documentation](https://htseq.readthedocs.io/en/master/htseqcount.html) for more details.

 The remaining optional parmeters allow you to tune and alter the functionality of all tools used by PROseq_etal. The top of this set includes three parameters that are probably best to check before running the pipeline. See config comments and linked documentation for details. The remaining are purely optional but can allow you to modify default settings of any tool used. **You never have to set parameters specifying output files or number of threads to be used**; PROseq_etal will handle these automatically.