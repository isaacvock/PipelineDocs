## Rule for counting mutations in a bam file

This is the big one, the custom scripts for counting mutations in a bam file. This is also a place where previous Simon lab post-doc Martin Machyna deserves a ton of credit for significantly improving the runtime and RAM efficiency of this step of the pipeline. Let's see how it works.

A brief summary of the key points:

1) A bam file and SNP calls are taken as input.
2) The bam file is fragmented into multiple smaller bam files, so that mutation counting can be parallelized using GNU parallel.
    - Sorting the bam file by query name was so as to facilitate fragmenting the bam files in this rule while keeping read pairs together in the same fragmented bam file.
3) Mutations of all types are tallied using a custom Python script. This script streams the bam file in one read at a time, and writes the mutation quantification to a csv file one read pair at a time. Thus, very little RAM is used in this step. 
    - Sorting the bam file by query name also is necessary to ensure proper mutation counting by this script. Mutation tallying is overwritten once a new read pair is arrived at, so if reads are not next to their pairs, the mutation counts for a read will be overwritten before it gets written to the final csv file.
    - In addition, it was important to filter out reads without their proper pair for the same reason: the mutation counting script will discard information for lone reads.
4) The `minqual` parameter in the config specifies the minimum quality score for a base call to contribute to the mutation tally. A nucleotide must also be a currently hardcoded minimum distance of 5 nts from the ends of the read to contribute to the mutation tally. The latter is due to the fact that base call quality gets worse at the ends of a read, but the current hardcoding is admittedly suboptimal and will be dealt with in later releases of bam2bakR.
5) Dovetailing (when the two reads of a pair include shared sequences) is resolved by giving precedence to higher confidence base calls in the dovetailed regions. If two base calls disagree and have the same quality score, the first read in the pair's base call is given precedence by default. This is an arbitrary decision akin to just randomly choosing between the discordant base calls.


### Snakemake rule (in workflow/rules/bam2bakR.smk)

``` python
# Count mutations
rule cnt_muts:
    input:
        "results/sf_reads/{sample}.s.bam",
        "results/snps/snp.txt",
    params:
        format=FORMAT,
        minqual=config["minqual"],
        mut_tracks=config["mut_tracks"],
        strand=config["strandedness"],
        shellscript=workflow.source_path("../scripts/mut_call.sh"),
        pythonscript=workflow.source_path("../scripts/mut_call.py"),
        awkscript=workflow.source_path("../scripts/fragment_sam.awk"),
        mutpos=False,
    output:
        "results/counts/{sample}_counts.csv.gz",
        temp("results/counts/{sample}_check.txt"),
    log:
        "logs/cnt_muts/{sample}.log",
    threads: 32
    conda:
        "../envs/full.yaml"
    shell:
        """
        chmod +x {params.shellscript}
        chmod +x {params.pythonscript}
        chmod +x {params.awkscript}
        {params.shellscript} {threads} {wildcards.sample} {input} {output} {params.minqual} {params.mut_tracks} {params.format} {params.strand} {params.pythonscript} {params.awkscript} {params.mutpos} 1> {log} 2>&1
        """
```

Ooh buddy that's a lot of custom scripts. First off, let's note the 5 non-script related parameters:

1) format: This will be either "SE" (for single-end) or "PE" (for paired-end).
2) minqual: This is an integer specified in the config file specifying the minimum base call quality score required for apparent mutations at a particular site to be called. 
3) mut_tracks: Not my best naming scheme, but refers to the types of mutations that will be included in the final cB, and for which color-coded coverage tracks will be made. This is specified in the config.
4) strand: "F" or "R", specified in the config. "F" means "forward stranded", i.e. the first (or only) read in a pair represents the original sequence of the RNA. "R" means "reverse stranded", i.e. the first (or only) read in a pair represents the reverse complemented sequence of the RNA from which it was derived.
5) mutpos: Eventually, this pipeline will support creating a table that tracks the mutational content of individual nucleotides. In fact, it is actually already implemented in the mutation counting script I will show later, as well as in the overhaul of this pipeline that currently is a separate repo ([fastq2EZbakR](https://github.com/isaacvock/fastq2EZbakR)). The vision is for bam2bakR (or fastq2EZbakR) to become a one-stop shop for processing all RNA-seq extensions which encode useful information as mutations in sequencing reads. Who knows how realistic that vision is...

Currently this rule caps the number of provided cores at 32. From some recent benchmarking, I have found that this rule makes highly effective use of even the max 32 cores. Thus, I will probably increase this cap in the near future. 

### Custom scripts called by cnt_muts

Before getting into specifics, let's layout a general view of what is happening in this rule:

1) The provided bam file is split up into multiple smaller bam files. This is done via a combination of a custom shell script and a custom AWK script. The reason for this is so that [GNU parallel](https://www.gnu.org/software/parallel/) can be used to count mutations in many of these smaller bam files in parallel.
2) Mutations are counted in the fragmented bam files, producing a table quantifying the number of occurences of every kind of mutation in each read. This is done with a custom Python script written originally by Martin Machyna and edited a bit by me (Isaac).
3) The tables generated by this step are combined into a single final table which is the desired output of this rule.

There are three custom scripts to discuss here. First is the shell script that gets initially called (as with other entries in the architecture page, some edits have been made to facilitate proper code highlighting and guide discussion):

``` bash
#!/bin/bash

# STEP 1:
# Source the paths and variables:
cpus=$1
sample=$2
input=$3
input_snp=$4
output=$5
output2=$6
minqual=$7
mut_tracks=$8
format=$9
strand=${10}
mutcnt=${11}
awkscript=${12}
mutpos=${13}

# Exit immediately if any command returns a non-zero status
set -e

# STEP 2:
# Create results/counts/
touch "$output2"

# Infer a decent fragment size
num_reads=$(samtools view -@ "$cpus" -c "$input")

fragment_size=$(echo "scale=0; $num_reads/$cpus" | bc)


# Make sure fragment size isn't 0, though if it were, you probably have bigger problems on your hand...
(( fragment_size++ ))

# STEP 3:
# Spread the work load so all cpus are working at all times
    declare $(samtools view -@ "$cpus"  "$input" 
    			| wc -l 
    			| awk -v fragment_size="$fragment_size" 
                      -v cpus="$cpus" '{
                                        nFrag = int($1 / fragment_size)                         # Calculate to how many fragments must the .bam be split
                                        if ($1 % fragment_size != 0) {nFrag = nFrag + 1}        # If there are some left over reads then increase number of fractions

                                        nBatchRun = int(nFrag / cpus)                           # In how many batches will the fragments be processed
                                        if (nFrag % cpus != 0) { nBatchRun = nBatchRun + 1}     # If there are some leftover fragments that will not fill up all cpus, then incerease number of batches

                                        newFragmentNumber = nBatchRun * cpus

                                        newFragmentSize = int($1 / newFragmentNumber)
                                        if ($1 % newFragmentNumber != 0) {newFragmentSize = newFragmentSize + 2}    # if there are still some reads leftover, then increase fragment size by 1 read-pair

                                        print "newFragmentNumber="newFragmentNumber
                                        print "newFragmentSize="newFragmentSize
                                    }')

    echo "* The number of fragments for sample $sample is $newFragmentNumber"
    echo "* The fragment size is set to $newFragmentSize alignments per fragment"



    for f in $(seq 1 $newFragmentNumber); do
        samtools view -H "$input" > ./results/counts/"$f"_"$sample".sam
    done &&


samtools view "$input" 
    | awk 
        -v fragment_size="$newFragmentSize" 
        -v sample="$sample" 
        -f "$awkscript" &&

for f in $(seq 1 $newFragmentNumber); do
    samtools view -@ "$cpus" -o ./results/counts/"$f"_"$sample"_frag.bam ./results/counts/"$f"_"$sample".sam
    rm ./results/counts/"$f"_"$sample".sam
done &&

echo "* Aligned .sam file fragmented for sample $sample"


# STEP 4:
# Call mutations
    parallel -j $cpus "python $mutcnt -b {1} \
                                              --mutType $mut_tracks \
                                              --minQual $minqual \
											  --SNPs "./results/snps/snp.txt" \
                                              --strandedness $strand \
                                              $( if [ "$mutpos" = "True" ]; then echo '--mutPos '; fi ) \
                                              --reads $format" ::: ./results/counts/*_"$sample"_frag.bam \


    echo "** Mutations called for sample $sample"


# STEP 5:
# Combine output from fragments into single file
    # 1) _count.csv files
    awk 'FNR > 1 || NR == 1' ./results/counts/*_"$sample"_frag_counts.csv 
        | pigz -p $cpus > "$output"

    rm ./results/counts/*_"$sample"_frag_counts.csv



    echo "** Results fragments merged into final files"



	rm -f ./results/counts/*_"$sample"_frag.bam

	echo '* Cleaning up fragmented .bam files'

```

There are effectively 5 steps in this script:

1) Parse the arguments provided by the Snakemake rule
2) Some preprocessing. Most noteable is the automatic setting of a "fragment_size". As mentioned above, the input bam file will be split up into multiple small bam files. There is no reason for the number of fragments to be greater than the number of provided cpus. Thus, a fragment size (number of reads per fragmented bam file) is calculated accordingly.
3) The input bam file is fragmented. This is accomplished with an AWK script written by Martin Machyna copied below. Read pairs are kept together, which is crucial for the mutation counting script. This is why sort_filter sorts the reads by query name. Doing so ensures that read pairs are right next to eachother in the bam file:

``` bash
#! /usr/bin/awk -f 
# Sep file into fragments 
# pass variables: fragment_size [-v fragment_size="300000"]
#                 sample [-v sample="my_sample_name"]

BEGIN { i = 1 }

# Remember read name for second last read
NR == (fragment_size * i - 1 ) { x = $1 }

# Keep adding reads to file until fragment size limit is reached
NR < (fragment_size * i ) { print >> "./results/counts/"i"_"sample".sam" }

# If the last read is a pair to previous one, add it to the same file. Othewise, start new fragment file.
NR == (fragment_size * i ) { if ($1 == x) 
                                { 
                                    print >> "./results/counts/"i"_"sample".sam"
                                    i++
                                } 
                            else 
                                {           
                                    i++
                                    print >> "./results/counts/"i"_"sample".sam"
                                }
                        
                            }

```

4) Mutations are counted using a custom Python script originally written by Martin Machyna and modified by me (Isaac). What is nice about this script is that the bam file is streamed. That means that the bam file is not loaded into memory all at once. Rather, it is read in via the pysam library one read at a time. Information about the mutation counts in that read are recorded and held in memory until the read's pair is encountered and its mutations counted. Once this happens, the read pair data is written to a csv file that grows one row at a time as the script runs, and the mutation count data is overwritten. Therefore, this script uses very little RAM. I will devote a whole section to the Python script, so its code will be reproduced later. The last thing to note is that GNU parallel is used to count mutations in multiple bam file fragments in parallel.
5) Finally, the mutation count tables are merged into a single table.


### The mutation counting script

A bulk of the interesting decisions are made in one Python script. I'll reproduce the full thing at the bottom of this page, but first, I would like to walk through it section by section, highlighting its key aspects and decisions.

Section 1: The setup (some details irrelevant to the current discussion are removed)

``` python
import pysam
import csv
import argparse
import datetime

# Parse commandline arguments
parser = argparse.ArgumentParser(description='This is python implementation of TimeLapse mutation calling')
requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument('-b', '--bam', type=str, required=True, metavar = 'in_file.bam',
                    help='Bam file to process')
parser.add_argument('--reads', default='PE', type=str, choices=['PE', 'SE'],
                    help='Type of mutation to record (default: PE)')
parser.add_argument('--minDist', default=5, type=int, metavar = '<int>',
                    help='Base distance from read-end filter (default: 5)')
parser.add_argument('--minQual', default=40, type=int, metavar = '<int>',
                    help='Base minimal quality filter (default: 40)')
parser.add_argument('--SNPs', help='Path to snp.txt file')
parser.add_argument('--strandedness', default='F', type=str, choices=['F', 'R'],
                    help='Is first read forward or reverse orientation? F = forward is the default.')
args = parser.parse_args()

# name without .bam suffix
inputName = args.bam.split('.bam')[0]           

if args.strandedness == 'F':
    strand_check = True
else:
    strand_check = False

# Initialize variables
firstReadName = ''

# Mutation type dictionary
muts = {'TA': 0, 'CA': 0, 'GA': 0, 'NA': 0, 'AT': 0, 'CT': 0, 'GT': 0, 'NT': 0, 'AC': 0, 'TC': 0, 'GC': 0, 'NC': 0, 'AG': 0, 'TG': 0, 'CG': 0, 'NG': 0, 'AN': 0, 'TN': 0, 'CN': 0, 'GN': 0, 'NN': 0}

# DNA code for comp and revcomp transformation
DNAcode={'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'N': 'N', 'a': 't', 'c': 'g', 't': 'a', 'g': 'c', 'n': 'n'}  

# Final table header
header = ['qname', 'nA', 'nC', 'nT', 'nG', 'rname', 'FR', 'sj', 'TA', 'CA', 'GA', 'NA', 'AT', 'CT', 'GT', 'NT', 'AC', 'TC', 'GC', 'NC', 'AG', 'TG', 'CG', 'NG', 'AN', 'TN', 'CN', 'GN', 'NN']


# Information to be extracted from a read
r_info = [''] + 4*[0] + 3*['']

# 
dovetail = []
MDstore = {}


# Load SNPs for filtering
snp = {}
snpFile = open(args.SNPs, 'r')
for line in snpFile:
    line = line.strip().split(':')
    snp[line[2] + ':' + line[3]] = line[0] + ':' + line[1]


#  Set .csv file for writing (simulating _counts.rds file)
myfile = open(inputName + '_counts.csv', 'w', newline='')
wr = csv.writer(myfile)
wr.writerow(header)


# Set .bam file for reading
samfile = pysam.AlignmentFile(args.bam, 'rb')
```

At the top of the script, arguments passed to the script are parsed, and a bunch of import variables are initialized. The arguments are:

1) `bam`: Path to the input bam file
2) `reads`: Bit of a weird name choice, but string denoting whether or not the reads are paired-end or single-end. This is important in the mutation counting portion, where a decision has to be made about whether to wait for the next read in the pair before writing to the final table or not.
3) `minDist`: This is the minimum distance that a nucleotide must be from the end of a read for a mutation to be counted. Base calling quality gets worse at the ends, so this filter was implemented. That being said, it probably shouldn't be a hard-coded magic number, and should rather be tunable in the config.
4) `minQual`: This is the minimum quality score for a base call to contribute to the mutation counting tally. This is set in the config
5) `SNPs`: This is the path to the snp.txt file, with information about the location of SNPs. This ensures that any tallied mutations are not just results of genomic differences between the reference and your cells/organism.
6) `strandedness`: This is a string representing whether or not the reads are forward or reverse stranded, and is set in the config.

Next, a number of variables are initialized, and relevant files created/opened. If you take your time to read the comments, it should be clear the roles of most of these variables/files. The `r_info`, `dovetail`, and `MDstore` objects might be a bit more cryptic, but their purpose will become clear as we move to the crux of the code, the actual mutation counting.

In the rest of the script, each read in the loaded bam file fragment is looped over and information about the base and mutation content of the read is tracked.


``` python
for r in samfile:

    # Initialize + acquire info: First read only
    if firstReadName != r.query_name:
        muts={'TA': 0, 'CA': 0, 'GA': 0, 'NA': 0, 'AT': 0, 'CT': 0, 'GT': 0, 'NT': 0, 'AC': 0, 'TC': 0, 'GC': 0, 'NC': 0, 'AG': 0, 'TG': 0, 'CG': 0, 'NG': 0, 'AN': 0, 'TN': 0, 'CN': 0, 'GN': 0, 'NN': 0}
        r_info = [''] + 4*[0] + 3*['']
        dovetail = []
        MDstore = {}
        gmutloc = []
        tp = []

        r_info[0] = r.query_name            # Read name
        r_info[5] = r.reference_name        # Chromosome name


    # Gather alignmet information + Resolve dovetailing: Both reads
    if ('I' not in r.cigarstring) and ('D' not in r.cigarstring):       # Any read without insertions/deletions

        r_info[7] = str( r_info[7] == 'TRUE' or ('N' in r.cigarstring) ).upper()     # sj: splice junction

        if (r.is_paired and (r.is_read1 == (r.is_reverse == strand_check))) or (not r.is_paired and (r.is_reverse == strand_check)):        # If read is first_in_pair and on reverse strand -or- second_in_pair and on forward strand then make sequence complement
            r_info[6] = 'R'      # FR: forward or reverse read orientation
            MD = [[x[1], DNAcode[x[2]], min(x[0] - r.query_alignment_start + 1, r.query_alignment_length - (x[0] - r.query_alignment_start))] for x in r.get_aligned_pairs(matches_only = True, with_seq=True)]
            # Parse MD and Cigar strings, remove values that are softclipped
            # MD = [[gen_position, ref_base, base_readEnd_distance]]

            temp_qual = r.query_qualities
            r.query_sequence = ''.join([DNAcode[x] for x in r.query_sequence])
            r.query_qualities = temp_qual
        else:
            r_info[6] = 'F'
            MD = [[x[1], x[2], min(x[0] - r.query_alignment_start + 1, r.query_alignment_length - (x[0] - r.query_alignment_start))] for x in r.get_aligned_pairs(matches_only = True, with_seq=True)]



        if firstReadName != r.query_name:       # First read
            MDstore = {z[0][0]: [z[0][1], z[1], z[2], z[0][2]] for z in zip(MD, r.query_alignment_sequence, r.query_alignment_qualities)}
            # store informatinon in dictionary of lists: {gen_position: [ref_base, read_base, qual, base_readEnd_distance]}
        else:                                   # Second read
            dovetail = list(set(MDstore.keys()) & set([x[0] for x in MD]))   # Identify genomic positions that are covered by both first and second read


            if len(dovetail) == 0:      # No dovetailing
                MDstore.update({z[0][0]: [z[0][1], z[1], z[2], z[0][2]] for z in zip(MD, r.query_alignment_sequence, r.query_alignment_qualities)})

            else:                       # Dovetailing
                MD = {z[0][0]: [z[0][1], z[1], z[2], z[0][2]] for z in zip(MD, r.query_alignment_sequence, r.query_alignment_qualities)}

                # MDstore.update({ pos:data for pos, data in MD.items() if pos in dovetail and MDstore[pos][2] < data[2] })   # Replace dovetail positions if better quality
                MDstore.update({ pos:data for pos, data in MD.items() if pos in dovetail and ((MDstore[pos][2] < data[2] and MDstore[pos][0].islower() and data[0].islower()) or (MDstore[pos][2] < data[2] and MDstore[pos][0].isupper() and data[0].isupper()) or (MDstore[pos][2] < data[2] and MDstore[pos][0].islower() and data[0].isupper() and MDstore[pos][2] + 33 < args.minQual) or (data[0].islower() and MDstore[pos][0].isupper() and data[2] + 33 > args.minQual)) })
                # This is a hack to simulate TimeLapse.R behaviour, but does not necessarily mean that it is a correct dovetail mutations handling
                # For dovetail bases: 1) If there is no mutation in 1st and in 2nd read => replace with higher quality 2nd read
                #                     2) If there is mutation in 1st and in 2nd read => replace with higher quality 2nd read
                #                     3) If there is mutation in 1nd but not in 2st read => replace only if 1st read quality is less than threshold and less than 2nd read
                #                     4) If there is mutation in 2nd read but not in 1st read => replace even with lower quality 2nd read as long as 2nd read quality is higher than threshold


                MDstore.update({ pos:data for pos, data in MD.items() if pos not in dovetail })                               # Append non dovetail sites



    # Collect data: Second read only or if in SE mode
    if (args.reads == 'SE' or firstReadName == r.query_name) and len(MDstore) > 0:
        refseq = [x[0].upper() for x in MDstore.values() if x[2] + 33 > args.minQual]  # Get reference sequence for readpair keeping only bases with given qaulity (Note: I think this should be also filtered for closeness to read end and presence of SNPs)
        # Count bases in reference sequence (soft clipped, dovetail-free)
        r_info[1] = refseq.count('A')       # nA
        r_info[2] = refseq.count('C')       # nC
        r_info[3] = refseq.count('T')       # nT
        r_info[4] = refseq.count('G')       # nG


        # Loop through every base of alignment and find mutations
        for pos, b in MDstore.items():

            # Find mutations marked as lowercase letters; apply quality filter; apply distance to read end filter; position is not a SNP
            if b[0].islower() and (b[2] + 33 > args.minQual) and (b[3] > args.minDist) and (r.reference_name + ':' + str(pos + 1) not in snp):   
                # Increment the mutation counter for current readpair
                muts[b[0].upper() + b[1]] += 1                                            




        # Write read info into _counts.csv
        r_info.extend( list(muts.values()) )
        wr.writerow(r_info)


    # Save read name for next iteration
    firstReadName = r.query_name

print('end: ' + str(datetime.datetime.now()))


##### Close files ######
myfile.close()
```


The full mutation counting script (unredacted), looks like:

``` python
'''
    File name: mut_call.py
    Author: Martin Machyna
    Email: machyna@gmail.com
    Orcid: 0000-0002-3624-3472
    Wikidata: Q55137016
    Date created: 8/20/2021
    Date last modified: 8/30/2021
    Version: 1.0.0
    License: GPLv3
    Python Version: Python 2.7.18, 3.8.7+
    Packages: pysam 0.16.0.1
    Description: Script for detecting mutations from sequencing data .bam file
'''

import pysam
import csv
import argparse
import datetime

# Parse commandline arguments
parser = argparse.ArgumentParser(description='This is python implementation of TimeLapse mutation calling')
requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument('-b', '--bam', type=str, required=True, metavar = 'in_file.bam',
                    help='Bam file to process')
parser.add_argument('--mutType', default='TC', type=str,
                    help='Type of mutation to record (default: TC)')
parser.add_argument('--reads', default='PE', type=str, choices=['PE', 'SE'],
                    help='Type of mutation to record (default: PE)')
parser.add_argument('--minDist', default=5, type=int, metavar = '<int>',
                    help='Base distance from read-end filter (default: 5)')
parser.add_argument('--minQual', default=40, type=int, metavar = '<int>',
                    help='Base minimal quality filter (default: 40)')
parser.add_argument('--tracks', action='store_true', # Automatically stored as default = FALSE
                    help='Generate files necessary for creating browser tracks')
parser.add_argument('--mutsRDS', action='store_true',
                    help='Generate _muts.rds like file with mutation frequency output')
parser.add_argument('--mutPos', action='store_true',
                    help='Generate _cU.rds like file and mutation position bedGraph files')
parser.add_argument('--SNPs', help='Path to snp.txt file')
parser.add_argument('--strandedness', default='F', type=str, choices=['F', 'R'],
                    help='Is first read forward or reverse orientation? F = forward is the default.')
args = parser.parse_args()

args.mutType = args.mutType.split(',')
args.base = [x[0] for x in args.mutType]        # base nucleotide: TC => T
inputName = args.bam.split('.bam')[0]           # name without .bam suffix

if args.strandedness == 'F':
    strand_check = True
else:
    strand_check = False



# Initialize variables
freq = {}               # dictionary for _muts.csv file with structure -> key='chrom:pos'  valuses=[trials, muts]
cU = {}
firstReadName = ''
muts = {'TA': 0, 'CA': 0, 'GA': 0, 'NA': 0, 'AT': 0, 'CT': 0, 'GT': 0, 'NT': 0, 'AC': 0, 'TC': 0, 'GC': 0, 'NC': 0, 'AG': 0, 'TG': 0, 'CG': 0, 'NG': 0, 'AN': 0, 'TN': 0, 'CN': 0, 'GN': 0, 'NN': 0}
DNAcode={'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'N': 'N', 'a': 't', 'c': 'g', 't': 'a', 'g': 'c', 'n': 'n'}  # DNA code for comp and revcomp transformation
header = ['qname', 'nA', 'nC', 'nT', 'nG', 'rname', 'FR', 'sj', 'TA', 'CA', 'GA', 'NA', 'AT', 'CT', 'GT', 'NT', 'AC', 'TC', 'GC', 'NC', 'AG', 'TG', 'CG', 'NG', 'AN', 'TN', 'CN', 'GN', 'NN']

# For counting mutations at individual positions
if args.mutPos:
    header.extend(['gmutloc', 'tp'])


r_info = [''] + 4*[0] + 3*['']
dovetail = []
MDstore = {}


# Load SNPs for filtering
snp = {}
snpFile = open(args.SNPs, 'r')
for line in snpFile:
    line = line.strip().split(':')
    snp[line[2] + ':' + line[3]] = line[0] + ':' + line[1]

# Create files for tracks
if args.tracks:
    # Create output files names
    fileName = []
    for mt in args.mutType:
        for i in range(0,6):
            fileName.append('_'.join([inputName, mt, str(i), 'reads.txt']))

    # Open all files for writing
    fs = []
    for f in fileName:
        fs.append(open(f, 'w'))



#  Set .csv file for writing (simulating _counts.rds file)
myfile = open(inputName + '_counts.csv', 'w', newline='')
wr = csv.writer(myfile)
wr.writerow(header)

# NO LONGER WRITING TO CU AS I ITERATE OVER BAM FILE
# TO SAVE ON DISK SPACE.
# # Initialize cU file
# if args.mutPos:
#     wr.writerow(['rname', 'gloc', 'GF', 'XF', 'ai', 'tp', 'trials', 'n'])
#
#     mycU = open(inputName + '_cU.csv', 'w', newline='')
#     cUwr = csv.writer(mycU)
#     cUwr.writerow(['rname', 'gloc', 'GF', 'XF', 'ai', 'tp', 'trials', 'n'])



# Set .bam file for reading
samfile = pysam.AlignmentFile(args.bam, 'rb')

print('Start: ' + str(datetime.datetime.now()))
for r in samfile:

    # Initialize + acquire info: First read only
    if firstReadName != r.query_name:
        muts={'TA': 0, 'CA': 0, 'GA': 0, 'NA': 0, 'AT': 0, 'CT': 0, 'GT': 0, 'NT': 0, 'AC': 0, 'TC': 0, 'GC': 0, 'NC': 0, 'AG': 0, 'TG': 0, 'CG': 0, 'NG': 0, 'AN': 0, 'TN': 0, 'CN': 0, 'GN': 0, 'NN': 0}
        r_info = [''] + 4*[0] + 3*['']
        dovetail = []
        MDstore = {}
        gmutloc = []
        tp = []

        r_info[0] = r.query_name            # Read name
        r_info[5] = r.reference_name        # Chromosome name


    # Gather alignmet information + Resolve dovetailing: Both reads
    if ('I' not in r.cigarstring) and ('D' not in r.cigarstring):       # Any read without insertions/deletions

        r_info[7] = str( r_info[7] == 'TRUE' or ('N' in r.cigarstring) ).upper()     # sj: splice junction

        if (r.is_paired and (r.is_read1 == (r.is_reverse == strand_check))) or (not r.is_paired and (r.is_reverse == strand_check)):        # If read is first_in_pair and on reverse strand -or- second_in_pair and on forward strand then make sequence complement
            r_info[6] = 'R'      # FR: forward or reverse read orientation
            MD = [[x[1], DNAcode[x[2]], min(x[0] - r.query_alignment_start + 1, r.query_alignment_length - (x[0] - r.query_alignment_start))] for x in r.get_aligned_pairs(matches_only = True, with_seq=True)]
            # Parse MD and Cigar strings, remove values that are softclipped
            # MD = [[gen_position, ref_base, base_readEnd_distance]]

            temp_qual = r.query_qualities
            r.query_sequence = ''.join([DNAcode[x] for x in r.query_sequence])
            r.query_qualities = temp_qual
        else:
            r_info[6] = 'F'
            MD = [[x[1], x[2], min(x[0] - r.query_alignment_start + 1, r.query_alignment_length - (x[0] - r.query_alignment_start))] for x in r.get_aligned_pairs(matches_only = True, with_seq=True)]



        if firstReadName != r.query_name:       # First read
            MDstore = {z[0][0]: [z[0][1], z[1], z[2], z[0][2]] for z in zip(MD, r.query_alignment_sequence, r.query_alignment_qualities)}
            # store informatinon in dictionary of lists: {gen_position: [ref_base, read_base, qual, base_readEnd_distance]}
        else:                                   # Second read
            dovetail = list(set(MDstore.keys()) & set([x[0] for x in MD]))   # Identify genomic positions that are covered by both first and second read


            if len(dovetail) == 0:      # No dovetailing
                MDstore.update({z[0][0]: [z[0][1], z[1], z[2], z[0][2]] for z in zip(MD, r.query_alignment_sequence, r.query_alignment_qualities)})

            else:                       # Dovetailing
                MD = {z[0][0]: [z[0][1], z[1], z[2], z[0][2]] for z in zip(MD, r.query_alignment_sequence, r.query_alignment_qualities)}

                # MDstore.update({ pos:data for pos, data in MD.items() if pos in dovetail and MDstore[pos][2] < data[2] })   # Replace dovetail positions if better quality
                MDstore.update({ pos:data for pos, data in MD.items() if pos in dovetail and ((MDstore[pos][2] < data[2] and MDstore[pos][0].islower() and data[0].islower()) or (MDstore[pos][2] < data[2] and MDstore[pos][0].isupper() and data[0].isupper()) or (MDstore[pos][2] < data[2] and MDstore[pos][0].islower() and data[0].isupper() and MDstore[pos][2] + 33 < args.minQual) or (data[0].islower() and MDstore[pos][0].isupper() and data[2] + 33 > args.minQual)) })
                # This is a hack to simulate TimeLapse.R behaviour, but does not necessarily mean that it is a correct dovetail mutations handling
                # For dovetail bases: 1) If there is no mutation in 1st and in 2nd read => replace with higher quality 2nd read
                #                     2) If there is mutation in 1st and in 2nd read => replace with higher quality 2nd read
                #                     3) If there is mutation in 1nd but not in 2st read => replace only if 1st read quality is less than threshold and less than 2nd read
                #                     4) If there is mutation in 2nd read but not in 1st read => replace even with lower quality 2nd read as long as 2nd read quality is higher than threshold


                MDstore.update({ pos:data for pos, data in MD.items() if pos not in dovetail })                               # Append non dovetail sites



    # Collect data: Second read only or if in SE mode
    if (args.reads == 'SE' or firstReadName == r.query_name) and len(MDstore) > 0:
        refseq = [x[0].upper() for x in MDstore.values() if x[2] + 33 > args.minQual]  # Get reference sequence for readpair keeping only bases with given qaulity (Note: I think this should be also filtered for closeness to read end and presence of SNPs)
        # Count bases in reference sequence (soft clipped, dovetail-free)
        r_info[1] = refseq.count('A')       # nA
        r_info[2] = refseq.count('C')       # nC
        r_info[3] = refseq.count('T')       # nT
        r_info[4] = refseq.count('G')       # nG


        # Loop through every base of alignment and find mutations
        for pos, b in MDstore.items():

            # _muts.rds data
            if args.mutsRDS:
                if (b[0].upper() in args.base):
                    if (r.reference_name + ':' + str(pos)) not in freq:
                        freq[r.reference_name + ':' + str(pos)] = [1, 0]                   # Initialize counter for position
                    else:
                        freq[r.reference_name + ':' + str(pos)][0] += 1                    # Increment read coverage counter for given genomic position

                if b[0].islower() and ((b[0].upper() + b[1]) in args.mutType):
                    freq[r.reference_name + ':' + str(pos)][1] += 1                        # Increment mutation counter for given genomic position

            # cU.rds trial data
            if args.mutPos:
                if (b[0].upper() in args.base):
                    whichMut = [mut for mut in args.mutType if mut[0] == b[0].upper()]     # Find out which mutation types use this reference base e.g. T -> TC, TG, TA, TN
                    for mt in whichMut:
                        key = r.reference_name + ':' + str(pos) + ':' + mt
                        if key not in cU:
                            cU[key] = [1, 0]
                        else:
                            cU[key][0] += 1

            # _counts.rds data
            if b[0].islower() and (b[2] + 33 > args.minQual) and (b[3] > args.minDist) and (r.reference_name + ':' + str(pos + 1) not in snp):   # Find mutations marked as lowercase letters; apply quality filter; apply distance to read end filter; position is not a SNP
                muts[b[0].upper() + b[1]] += 1                                            # Increment the mutation counter for current readpair

                # mutPos bedGraph data + cU.rds n data
                if args.mutPos:
                    if (b[0].upper() + b[1]) in args.mutType:
                        key = r.reference_name + ':' + str(pos) + ':' + b[0].upper() + b[1]
                        cU[key][1] += 1

                        key = r.reference_name +  ':' + str(pos) + ':' + r_info[6] + ':' + b[0].upper() + b[1]

                        gmutloc.append(str(pos))            # Record position of muatation
                        tp.append(b[0].upper() + b[1])      # Record type of mutation




        # Write read info into _counts.csv
        r_info.extend( list(muts.values()) )
        if args.mutPos:
            r_info.extend( [ '|'.join(gmutloc), '|'.join(tp) ] )
        wr.writerow(r_info)

        # Was doing this to save on RAM
        # Realized this inflated disk space usage tremendously
        # Will save here for posterity's sake
        # # Write to cU file
        # if args.mutPos:
        #
        #     for position, counts in cU.items():
        #         row = position.split(':')
        #         row[1] = int(row[1]) + 1                        # adjust position because we are 0-based
        #         row.extend(counts)
        #         cUwr.writerow(row)
        # 
        #     cU = {}

        # Save read name to track files
        if args.tracks:
            for index, mut in enumerate(args.mutType):
                for c in range(0, (5 if muts[mut] > 5 else muts[mut]) + 1):
                    fs[c + index * 6].write(r.query_name + '\n')

    # Save read name for next iteration
    firstReadName = r.query_name

print('end: ' + str(datetime.datetime.now()))


##### Close files ######
myfile.close()

### saving cU.rds file
if args.mutPos:
    with open(inputName + '_cU.csv', 'w', newline='') as myfile:
        wr = csv.writer(myfile)
        wr.writerow(['rname', 'gloc', 'tp', 'trials', 'n'])
        for position, counts in cU.items():
            row = position.split(':')
            row[1] = int(row[1]) + 1                        # adjust position because we are 0-based
            row.extend(counts)
            wr.writerow(row)

    del cU
print('cU: ' + str(datetime.datetime.now()))



if args.tracks:
    for f in fs:
        f.close()
```

