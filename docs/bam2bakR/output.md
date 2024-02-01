## bam2bakR Output

All output files will be placed in a directory named `results` that will be created the first time you run bam2bakR. The output of greatest interest, the gzipped cB.csv file, will be in `results/cB/`. The cB file is a tidy table where each row corresponds to data for a set of reads with identical values for the following columns:

  * sample (sample ID)
  * GF: Gene ID for reads aligning anywhere in a gene (intron or exon)
  * XF: Gene ID for reads aligning to exclusively exonic regions of a gene
  * nT: Number of Ts in the reads (or As if the read comes from the reverse complement of the original RNA)
  * TC: Number of T-to-C mutations in the reads

The final column of the cB file is n, which represents the number of reads with the same value for sample, GF, XF, nT, and TC.

**Other bam2bakR v3.0.0+ output includes**

Processed bam files:

* Sorted and filtered bam/sam files are in `results/sf_reads/`
  - These are passed to the mutation counting and feature assignment scripts
  - Sorting and filtering is accomplished with a custom shell script.

Mutation counting output:

* `<sampleID>_counts.csv.gz` files are in `results/counts/`
  - Each row of this table corresponds to a single read or read pair
  - The columns represent counts of every mutation type and every type of nucleotide
  - These can be useful for tracking down problems in the mutation counting. See [FAQs](faqs.md) for details.
  - Mutation counting is accomplished with a custom python script called by a custom shell script.

Feature assignment with featureCounts

* Tables of exonic read counts for each annotated gene are in `results/featurecounts_exons/<sampleID>.featureCounts`.
* Table of detailed read assignment to exonic regions (in featureCounts' CORE format, a tsv file where each row provides information for where each read is assigned) are in `results/featurecounts_exons/<sampleID>.s.bam.featureCounts`.
* Table of number of reads supporting each exon-exon junction are in `results/featurecounts_exons/<sampleID>.jcounts`.
* Similar tables for assignment of reads to anywhere in a gene are in `results/featurecounts_genes`.


Merged feature assignment and mutation counting:

* Tables that have combined the exonic and gene feature assignment information with the mutation calling output are in `results/merge_feature_and_muts/<sampleID>_counts.csv.gz`. 
  - If a read was not assigned to a particular feature type (i.e., exon or gene), then it will have an NA in the relevant feature column (XF for exons and GF for genes).


Colored sequencing tracks:

* .tdf files that can be used to make the sequencing tracks colored by mutational content (described [here](../tracks.md)) are in `results/tracks`
  - Each sample has 12 .tdf files, named like `<sampleID>.TC.<#>.<strand>.tdf`, where `<#>` represents a number from 0 to 5 (number of T-to-C mutations in the reads used to make that file) and `<strand>` is either `pos` (plus strand) or `min` (minus strand).
  - Currently, if your library is reverse stranded (i.e., first read in a pair represents reverse complement of original RNA sequence), then the plus and minus strand tracks will be flipped. This does not change interpretation of the tracks, you just have to be aware of that when using an annotation to visually decide what reads are the 
  product of sense and antisense transcription.

Single nucleotide polymorphism (SNP) calls:

* SNP calls are in two formats (.txt and VCF)in the `results/snps/` directory.
  - If you did not have any -s4U control samples, then the .vcf file will not exist and the .txt file will be empty
  - These SNP calls are used to identify nucleotides which should be ignored for T-to-C mutation counting

Normalization:

* Scale factors calculated using edgeR's TMM strategy are located in `results/normalization/scale`
  - This is a simple tab-delimited text file with two "columns", one corresponding to the sample ID, and the other corresponding to the scale factor
  - These will be used to scale the heights of the sequencing tracks in `results/tracks/`.


**Additional output if providing fastq's as input**

Trimmed fastq files:

* Cutadapt output is in `results/fastq_cut/`

Bam files:

* Alignments provided by whatever aligner you chose are in `results/bams/`