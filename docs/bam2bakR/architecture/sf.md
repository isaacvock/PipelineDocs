## Rule for sorting and filtering a bam file

This is most likely the first rule that will execute when you run bam2bakR. It sorts and filters an input bam file using SAMtools and some custom AWK code. Let's see how this works.

### Snakemake rule (in workflow/rules/bam2bakR.smk)

``` python
# Filter out multi-mappers and sort reads
rule sort_filter:
    input:
        get_input_bams,
    output:
        "results/sf_reads/{sample}.s.bam",
        "results/sf_reads/{sample}_fixed_mate.bam",
        "results/sf_reads/{sample}.f.sam",
    log:
        "logs/sort_filter/{sample}.log",
    params:
        shellscript=workflow.source_path("../scripts/sort_filter.sh"),
    threads: 8
    conda:
        "../envs/full.yaml"
    shell:
        """
        chmod +x {params.shellscript}
        {params.shellscript} {threads} {wildcards.sample} {input} {output} {config[FORMAT]} 1> {log} 2>&1
        """
```

See [Tips and Tricks](../../pragmatism.md) for a description of all parts of a Snakemake rule. Pretty much every rule in bam2bakR looks similar to this. The input in this case is a function defined in `workflow/rules/common.smk` as:

``` python
def get_input_bams(wildcards):
    return config["samples"][wildcards.sample]
```

get_input_bams thus returns the specified path to the bam file provided by the user. The input to this function is a [Snakemake wildcard](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html), and the `sample` wildcard is the sample ID specified in the config file by the user.

A single conda environment is specified for all rules in bam2bakR. This is a bit sloppy and will likely change in the future, but it reduces the amount of time spent creating environments when you first launch the pipeline (though the runtime benefit is admittedly pretty trivial). In addition, this detail has little impact on the user experience.

This rule will use at most 8 threads. This is due to the fact that under the hood (see below), samtools will be called, and samtools rarely benefits from more than 8 threads for the tasks at hand.

This rule specifies the path to a script that will be called in the `shell` block. This script has to be made executable (that is the first line in the `shell` block) before calling it (second line). Most of the rules in bam2bakR call custom scripts like this. This leads to a somewhat annoying quirk in newer versions of Snakemake that means by default, Snakemake will always rerun these rules, because it considers the copying of the required script to be a "change" since the last time you ran the pipeline. See [here](https://github.com/snakemake/snakemake/issues/1694) for a further discussion of this. This is why in the documentation I suggest running Snakemake with the option `--rerun-triggers mtime`, which will only rerun rules if the input has been modified.

### Custom script called by Snakemake rule

This is one of many custom scripts used by bam2bakR. In addition, it is one of many legacy scripts in bam2bakR that I (Isaac Vock) did not originally write. Who wrote this script originally? Not sure, but it had to be one of: Matthew Simon, Jeremy Schofield, Martin Machyna, Lea Kiefer, or Joshua Zimmer. The versions of legacy scripts used in bam2bakR are most likely attributable mainly to Martin Machyna, who was a previous post-doc in the Simon lab who put an immense amount of effort into optimizing the pipeline. The major improvements he made will be called out in the various pipeline architecture links. The original form of these scripts can be found in a [bitbucket](https://bitbucket.org/mattsimon9/timelapse_pipeline/src/master/) housing the old, non-Snakemake TimeLapse pipeline.  

The script used in bam2bakR is reproduced below. Some carriage return escapes have been removed to get the code highlighting to cooperate. I have also added some comments to guide the discussion below:

``` bash
#!/bin/bash

# STEP 1:
# Source the paths and variables:
cpus=$1
sample=$2
input=$3
output=$4
output2=$5
output3=$6
format=$7

# STEP 2:
# Make sure bam file is query sorted rather than coordinate sorted
samtools sort -@ "$cpus" -n "$input" | samtools fixmate -@ "$cpus" - - | samtools view -@ "$cpus" -b - > "$output2"

# STEP 3:
	if [ "$format" = "NU" ]; then
        samtools view -@ "$cpus" -h "$output2" 
            | awk '$1 ~ /^@/ {print}
                    (($2 == 147 || $2 == 99) || ($2 == 83 || $2 == 163)) || (($2 == 355 || $2 == 403) || ($2 == 339 || $2 == 419))  {print}' > "$output3"
    elif [ "$format" = "PE" ]; then
	    samtools view -@ "$cpus" -q 2 -h "$output2" 
            | awk '$1 ~ /^@/ {print}
                  ($2 == 147 || $2 == 99) || ($2 == 83 || $2 == 163) {print}' > "$output3"
    elif [ "$format" = "SE" ]; then
        samtools view -@ "$cpus" -q 2 -h "$output2" 
            | awk '$1 ~ /^@/ {print}
                  ($2 == 0 || $2 == 16) {print}' > "$output3"
	else
		samtools view -@ "$cpus" -h "$sample"_fixed_mate.bam > "$output3"
	fi &&
    echo "* Reads filtered for sample $sample"
	


# STEP 4:
samtools sort -@ "$cpus" -n -O bam -o "$output" "$output3"


```

This script can be broken into four steps, as commented in the above code:

1) Parsing the input parameters
    - there are better ways to do this in newer versions of Snakemake, but it works and its how I originally figured out how to pass Snakemake rule parameters/input/output to a custom script  
    - cpus represents the number of threads provided
    - sample is the sample ID (also the sample wildcard in the Snakemake rule)
    - input and the three outputs are the input bam file and the three output sam/bam files produced by the script
    - format specifies whether or not the bam file is paired-end or single-end.

2) Fix mates
    - Technically this is only relevant for paired-end data.
    - Makes sure that each read in a pair has correct information about its mate (the other read in the pair).
    - `samtools fixmate` requires query sorted bam files as input, so this sorting is done prior to running `fixmate`.
    - This creates the `_fixed_mate.bam` file output

3) Filter reads
    - I have yet to memorize the various SAM flags and what they mean. I always go to [this lovely resource](https://broadinstitute.github.io/picard/explain-flags.html) to remind myself.
    - Let's ignore the NU format, which isn't really supported but is here for historical purposes.
    - Reads with a mapping score less than 2 are ignored and filtered out. Why? The cutoff is a bit arbitrary, but let's just call these low confidence alignments. They are more likely to have spurious, alignment artifact mutations.
    - If the reads are paired-end, then 147 and 99 represent reads with a proper mate and that have been mapped to a region near where their mate has been mapped.
    - If the reads are single-end, then the any read that is mapped will be kept. The "primary" alignment of multi-mapped reads are also kept.

4) Sort reads by coordinate
    - I am partaking on this exercise for two reasons. One, I want to provide detailed documentation of every decision made inside of this pipeline, so that curious users can easily figure out what is going on in any given rule and why those decisions were made. Two, I want to take an opportunity to scrutinize parts of legacy scripts that I haven't had reason to think about carefully yet.
    - This is a case where the second motivation for this exercise is relevant. Aren't the reads already query sorted? Perhaps this step is redundant, or perhaps mate fixing and filtering with AWK causes the ordering of reads to get mixed up. I need to look into this...


### Future improvements

1) It's unclear why AWK is being used to filter reads with certain SAM flags. I am pretty sure this can all be done with samtools, and I suspect it would be more efficient. Maybe I am wrong about the efficiency thing, but I will need to test that.
2) Need to look into removing the redundant sorting of reads (and whether or not it is actually redundant).
3) Creating and saving multiple intermediate bam files is just a waste of disk space. Let's not do that.


