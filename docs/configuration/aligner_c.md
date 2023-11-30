## Configuring THE_Aligner

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
`samples` is the list of sample IDs and paths to .bam files that you want to process. Delete the existing sample names and paths and add yours. The sample names in this example are `WT_1`, `WT_2`, `WT_ctl`, `KO_1`, `KO_2`, and `KO_ctl`. These are the sample names that will append many of the files output by THE_Aligner. The `:` is necessary to distinguish the sample name from what follows, the path to the relevant bam file. Note, the path can be absolute (e.g., ~/path/to/fastqs/) or relative to the directory that you deployed to (i.e., `workdir` in this example). In the example config, the paths specified are relative. Thus, in this example, the bam files are located in a directory called `samples` that is inside of a directory called `data` located in `workdir`. Your data can be wherever you want it to be, but it might be easiest if you put it in a `data` directory inside the THE_Aligner directory as in this example. 

As another example, imagine that the `data` directory was in the directory that contains `workdir`, and that there was no `samples` subdirectory inside of `data`. In that case, the relative paths would look something like this:

``` yaml
samples:
  WT_1: ../data/WT_1
  WT_2: ../data/WT_2
  WT_ctl: ../data/WT_ctl
  KO_1: ../data/KO_1
  KO_2: ../data/KO_2
  KO_ctl: ../data/KO_ctl
```
where `../` means navigate up one directory. 

The remaining required parameters are:

* `PE`: True if fastqs are paired-end. False if they are single-end.
* `genome`: Path to genome fasta file to align to.
* `annotation`: Path to annotation gtf file. Used to generate splice-aware alignment indices or to quantify with respect to if using a pseudoaligner.
* `transcriptome`: Path to transcriptome fastq file. Only relevant for pseudoaligners, and will get created automatically at the specified path if it does not already exist.
* `aligner`: Determines which aligner will be used (options are `"bwa-mem2"`, `"bowtie2"`, `"star"`, `"hisat2"`,`"kallisto"`, `"salmon"`).
* `indices`: Path to aligner indices. These will be created at this path if not already present.
* `strandedness`: Parameter specifying library strandedness. Options are "reverse", "yes", or "no". Naming convention comes from HTSeq, but the parameter is used in THE_Aligner by Salmon and Hisat2.
* `directionality`: Whether the reads face inwards, outwards, or in the same direction. Only used by Salmon.

 The remaining optional parmeters allow you to tune and alter the functionality of all tools used by THE_Aligner. The top of this set includes one parameter that is probably best to check before running the pipeline; the adapters to trim. See config comments and linked documentation for details. The remaining are purely optional but can allow you to modify default settings of any tool used. **NOTE: You never have to set parameters specifying output files or number of threads to be used**; THE_Aligner will handle these automatically.
