## The architecture of bam2bakR

Here I will provide extensive information as to the inner workings of every step of bam2bakR. Like every Snakemake workflow, bam2bakR can be decomposed into a set of "rules", i.e. steps that accept some input and generate some output. Each rule has its own page linked below, where the code relevant to that rule will be posted and described. 

### bam2bakR version 3.0.0 rules (bam file input)

[sort_filter](architecture/sf.md): Bam files provided by the user are sorted and filtered prior to mutation counting and feature assignment.
    - Input: provided bam files.
    - Output: 1 sam file and 2 different bam files. the sam file and one of the bam files are processing intermediates preserved for posterity's/debugging's sake. The bam file with suffix `.s.bam` is used in all downstream steps.

[call_snps](architecture/snp.md): SNPs are identified using bam files from -s4U control samples
    - Input: Genome FASTA, Genome index, sorted and filtered bam files from -s4U control samples.
    - Output: SNP calls in .txt and .vcf files.

[cnt_muts](architecture/mut.md): Mutations of every type are identified and quantified in processed bam files
    - Input: sorted and filtered bam file, SNP calls.
    - Output: table of the number of occurences of every type of mutation and reference nucleotide present in each read in the processed bam file.

[featurecounts_genes](architecture/gene.md): Reads are assigned to annotated genes with featureCounts
    - Input: Sorted and filtered bam file, reference annotation in GTF format.
    - Output: Various featureCounts output tables, most notably one in the CORE format, which provides detailed information of each read's assignment status.

[featurecounts_exons](architecture/exons.md): Reads are assigned to annotated exons with featureCounts
    - Input: Sorted and filtered bam file, reference annotation in GTF format.
    - Output: Various featureCounts output tables, most notably one in the CORE format, which provides detailed information of each read's assignment status.

[normalize](architecture/norm.md): Normalization scale factors are calculated using edgeR; these will be provided to the rule creating coverage tracks (maketdf).
    - Input: featureCounts summary of reads mapping to each gene (exonic regions only).
    - Output: text file with scale factors for each sample.

[merge_features_and_muts](architecture/merge.md): Merge the tables of read assignments from featureCounts with the mutation count table.
    - Input: featureCounts CORE format tables from gene and exon assignments, mutation calling output tables.
    - Output: one table with feature assignment and mutation calling information for each read/read pair.

[maketdf](architecture/tdf.md): Create .tdf coverage tracks color-coded by T-to-C mutational content of reads.
    - Input: processed bam file, mutation count table, normalization scale factors.
    - Output: .tdf files for making color-coded coverage tracks


[makecB](architecture/cB.md): Create cB.csv.gz file that is the main desired output of bam2bakR.
    - Input: Merged feature assignment + mutation counting tables for all samples
    - Output: cB.csv.gz