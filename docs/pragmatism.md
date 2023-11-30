# Tips and Snakemake Information

This page includes some suggestions for how to work with Snakemake pipelines. It also introduces some basic Snakemake concepts, such as rules, wildcards, and how the desired final output is specified, with examples of all of these from the PROseq_etal pipeline.

## Pipeline usage tips

* Save a copy of your config file along with the output of the pipeline. This can facilitate rerunning the pipeline on datasets you already ran it on once. It also provides a record of the parameters you used for each pipeline run.
* If the pipeline fails due to time out of a rule, or because you purposely cancelled its jobs prematurely, the directory in which you ran the pipeline can become "locked". You'll know when this happens because trying to rerun the pipeline will yield an error saying something along the lines of "the directory is locked". To unlock it, activate an environment with Snakemake installed and run `snakemake --unlock` inside the locked directory. 
* If something goes wrong, make sure to check both the .log and .out files. The former usually contains the output of running a particular tool (e.g., bwa-mem2), and this is often most useful for identifying the source of a problem. In some cases though, such logs of the output of a tool are not possible to capture, leaving the .out file to capture any relevant error messages.
* **FOR SIMON LABBERS OR YALE USERS** If running the pipeline within scratch60, parts of the conda environments created by the pipeline can get deleted in the 60 day time period. You can tell this is the case when log files suggest that the relevant software is not available (e.g., "samtools: command not found"). You can force the pipeline to recreate the conda environemnts by deleting the hidden `.snakemake/` directory contained in your working directory (the directory in which you ran the pipeline): `rm -r .snakemake`
* One plus of using Snakedeploy to deploy pipelines is that you don't have to manually update the pipeline. If you are tracking the main branch of the repository, then you will always be running the most up-to-date code on that branch. There can be a bit of a lag (a couple minutes at the most) between updates being made to a branch, and those changes being registered by a Snakedeploy deployed workflow. Be mindful of this if you are trying to run the pipeline with recently made changes.
* At any time, you can change which branch or tag (archived previous version) of the pipeline you are using by going into the `workflow/` directory created when you first deploy the workflow with Snakedeploy. Inside, you will find a lone, rather concise Snakefile. The meat of it will look like:

``` python
module PROseq_etal:
    snakefile:
        github("isaacvock/PROseq_etal", path="workflow/Snakefile", branch = "main")
    config:
        config
```

>You can change `branch = "main"` to whatever existing branch you please.

## Structure of the pipeline

Pipelines that I developed (bam2bakR, PROseq_etal, THE_Aligner, etc.) are structured [as recommended](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#integrated-package-management) by the Snakemake developers. A Snakefile is located in the `workflow/` directory. This Snakefile is fairly barebones and merely indicates what the expected final output is, and what files contain all of the rules to run. Each step of the pipeline is specified in .smk files in the `workflow/rules/` directory. These files are split up into sets of similar rules (e.g., all of the rules for running MACS2 in PROseq_etal are in one .smk file, called `macs2.smk`).

Some of the rules are wrappers from the [Snakemake wrapper repository](https://snakemake-wrappers.readthedocs.io/en/stable/). These wrappers load remote conda environment specifications stored in a common Github repo. Some of the rules are custom made and thus require custom conda environments that are specific to PROseq_etal. These conda environments are specified in .yaml files in the `workflow/envs/` directory. They specify the dependencies for sets of related rules, and what conda channes they can be downloaded from.

Finally, some steps use custom scripts that I wrote. These scripts are located in the `workflow/scripts` directory. Python scripts are easy to call in Snakemake rules, but a little bit of extra effort is required to call non-Python scripts. Such rules typically look like:

``` python
params:
    rscript=workflow.source_path("../scripts/calculate_PI.R")
threads: 1
shell:
    r"""
    chmod +x {params.rscript}
    {params.rscript} -p {input.pause} -g {input.gb} -a {input.gtf} -o {output.PI}
    """
```

This is an example from the step of the PROseq_etal pipeline that calculates pause indices. The script has to be specified using `workflow.source_path(...)`, which causes Snakemake to create a temporary copy of the script on your system. This causes problems in newer versions of Snakemake, as the act of copying these scripts is often logged as a change that causes the relevant step to get rerun even if its inputs have remained unchanged. In the `run_slurm.sh` script in the yale_profile repository discussed in the Yale deployment documentation, that is why `--rerun-triggers mtime` is included in the call to Snakemake; this makes it so that bona fide modification of inputs is the only thing that will trigger a rerun a rule. Finally, note that the script needs to be made executable with `chmod +x {params.rscript}`.


## Anatomy of a Snakemake rule

Each step of a Snakemake pipeline is referred to as a "rule". Here, I will show you some example rules from PROseq_etal and walk through the syntax of the rule and how it works. Most rules will have a basic structure identical to this example rule:

``` python
### Create fold enrichment track
rule macs2_enrichment:
    input: 
        treatment="results/macs2_callpeak/{treatment}_treat_pileup.bdg",
        control="results/macs2_callpeak/{treatment}_control_lambda.bdg",
    output:
        fe="results/macs2_enrichment/{treatment}_FE.bdg"
    params:
        extra=config["bdgcmp_FE_params"],
    log:
        "logs/macs2_enrichment/{treatment}.log"
    conda:
        "../envs/macs2.yaml"
    threads: 4
    shell:
        """
        macs2 bdgcmp \
            -t {input.treatment} \
            -c {input.control} \
            -o {output.fe} \
            -m FE {params.extra} 1> {log} 2>&1 
        """
```

This rule runs MACS2 bedGraph comparison, calculating the fold enrichment of an enriched ChIP-seq sample to the unenriched Input sample. All rules are a given a name, specified in the first line `rule <rule name>:`. 

The next three lines specify the file that the rule expects as input; this is known as the "input block". Each subsection within a rule that has a name followed by a colon is called a block. In this case, there are two expected input files: one located at `"results/macs2_callpeak/{treatment}_treat_pileup.bdg"` and another located at `"results/macs2_callpeak/{treatment}_control_lambda.bdg"`. File paths are relative to the directory from which the pipeline was run These are both bedGraph files generated by a separate rule. Commas need to be included at the end of each separate file specification. The comma at the end of the second input file specification is not strictly necessary, but can help avoid bugs if additional input is added at a later date. This rule will not run until its input is created, and the step will never run if its output is not required to generate the final desired output of the entire pipeline. More on how the "final desired output of the entire pipeline" is specified later.

The output block specifies the file (or in some cases multiple files) that this rule is expected to produce. One thing you may have noticed in both the input and output file names is the odd `{treatment}` in the middle of the file path. This is called a wildcard and is a unique and highly useful feature of Snakemake. They can also be the greatest source of annoyance when first learning Snakemake. Snakemake does some clever work to infer all of the different values for `{treatment}` to expect, and will run this rule separately for these different related sets of input. In this case, `{treatment}` refers to enrichment sample IDs. You can thus replace `{treatment}` with the sample IDs included as key names in the `controls` block of your config file. In the example config provided when you deploy the pipeline with Snakedeploy, this would be `WT_1`, `WT_2`, `KO_1`, and `KO_2`. So the use of wildcards in this rule tells Snakemake that this rule needs to be run separately (and if possible in parallel) for the following input:

* `"results/macs2_callpeak/WT_1_treat_pileup.bdg"` and `"results/macs2_callpeak/WT_1_control_lambda.bdg"`
* `"results/macs2_callpeak/WT_2_treat_pileup.bdg"` and `"results/macs2_callpeak/WT_2_control_lambda.bdg"`
* `"results/macs2_callpeak/KO_1_treat_pileup.bdg"` and `"results/macs2_callpeak/KO_1_control_lambda.bdg"`
* `"results/macs2_callpeak/KO_2_treat_pileup.bdg"` and `"results/macs2_callpeak/KO_2_control_lambda.bdg"`

and will generate the following output:

* `"results/macs2_enrichment/WT_1_FE.bdg"`
* `"results/macs2_enrichment/WT_2_FE.bdg"`
* `"results/macs2_enrichment/KO_1_FE.bdg"`
* `"results/macs2_enrichment/KO_2_FE.bdg"`

The params block specifies any additional parameters you want to specify for the tool(s) that you will be running in this rule. In this case, `config["bdgcmp_FE_params"]` means to use the value specified in the `bdgcmp_FE_params` parameter of the user's config file. In PROseq_etal, I make extensive use of this block to allow users to tune the performance of every single tool used by the pipeline.

The log block specifies the location of a file that is meant to capture the output of any tool run in this rule. This includes all messages output by the tool as it is running, and any errors that crop up.

The threads block specifies the maximum number of "threads" (basically cpus) that this rule will make use of. If you provide more cpus to this rule, only this many will be used, allowing the remaining cpus to be used for other jobs. In this case, MACS2 is not explicitly multithreaded, but is able to make use of multiple CPUs due to the BLAS library used by the numpy package.

The conda block specifies the location of a yaml file specifying a conda environment in which to run this step. If you use the `--use-conda` option when running Snakemake, this will cause all such conda environments to be automatically created the first time you run the pipeline. This helps ensure reproducibility of the pipeline.

Finally there is the shell block, which specifies the shell code you want to run. The shell block is not strictly necessary and can be placed with a scripts block to run Python scripts, or a run block to run raw Python code. The code inside this block does one thing: runs MACS2's bedGraph comparison tool. Information in each of the blocks can be used in the shell block using `{block_name.item_name}`, where `block_name` in this case can either be `input`, `output`, `params`, `log`, `conda` (though pretty much never this), and `{threads}`. If a block has multiple named items (e.g., `treatment=...`), the the name specified to the left of the equal sign can be supplied as the `item_name`. The last bit of the shell code (`1> {log} 2>&1`) says to write stdout and stderr to the specified log file, which captures all messages created by running this tool.

In PROseq_etal, I have made extensive use of the Snakemake wrapper repository, a collection of pre-specified Snakemake rules that you can easily insert into your own pipeline to make use of a number of popular bioinformatic tools. A rule making use of a wrapper looks like:

``` python
# Run fastqc on trimmed fastqs
rule fastqc:
    input:
        "results/trimmed/{sample}.{read}.fastq"
    output:
        html="results/fastqc/{sample}_r{read}.html",
        zip="results/fastqc/{sample}_r{read}_fastqc.zip"
    log:
        "logs/fastqc/{sample}_r{read}.log"
    params:
        extra = config["fastqc_params"]
    resources:
        mem_mb = 9000 
    threads: 4
    wrapper:
        "v2.2.1/bio/fastqc"
```

In this particular case, the main difference is that the shell block has been replaced by a wrapper block, which specifies the exact version of the wrapper you want to use. Implicitly, this also defines a conda environment to create upon the first run of the pipeline. As a fun aside, this rule includes a second unique features, which is a resources block that sets a max on the amount of RAM that can be requested for this step.

I mentioned that a rule will only be run if its output is "required to generate the final desired output of the entire pipeline". How is this determined though? Snakemake has to be provied a so-called "target rule", and the pipeline will run all rules necessary to generate the input for that target rule. The pipeline will finish when this target rule's output is generated. By default, the first rule in the main Snakefile is considered the target rule. Therefore, the first rule defined in the Snakefile located in the `workflow/` directory of the pipeline looks like (at least for the ChIP-seq option):

``` python
rule all:
    input: 
        expand("results/fastqc/{SID}_{read}.html", SID = SAMP_NAMES, read = READ_NAMES),
        expand("results/fastqc/{SID}_{read}_fastqc.zip", SID = SAMP_NAMES, read = READ_NAMES),
        expand("results/align/{SID}.bam", SID = SAMP_NAMES),
        expand("results/bigwig/{SID}.bw", SID = SAMP_NAMES),
        expand("results/macs2_diff_bw/{TID}_diff.bw", TID = TREATMENT_NAMES),
        expand("results/macs2_FE_bw/{TID}_FE.bw", TID = TREATMENT_NAMES),
        expand("results/homer_annotatePeaks/{SID}_annot.txt", SID = SAMP_NAMES),
        "results/homer_annotatePeaks/merged_annot.txt"
```

This rule is much simpler than the other example I showed. It only has an input block and nothing else! This is a common trick to essentially say "these are all of the files I want to generate". Therefore, Snakemake will run whatever rules are necessary to generate these files. It effectively starts from these desired inputs and works backwards, figuring out which steps need to be run to achieve this final goal. 

Another important feature of this rule is the use of Snakemake's `expand()` function. Expand will output a list of file paths, in this case replacing any `{...}` with the respective value specified in the comma separated statements after the file path. I mentioned that Snakemake has to cunningly decipher what the values of any wildcards in a given rule are. The trick to making this work for PROseq_etal relies on the use of `expand` in this target rule. The `SAMP_NAMES`, `READ_NAMES`, and `TREATMENT_NAMES` parameters come from input provided by the user (plus a bit of simple Python scripting), and are the key to specifying the value of all wildcards present in the pipeline. `SAMP_NAMES` for example is the list of sample IDs specified in the config's `samples` section. Therefore, `expand("results/bigwig/{SID}.bw", SID = SAMP_NAMES)` specifies that the following output must be generated (asssuming the example sample IDs of `WT_1`, `WT_2`, `WT_ctl`, `KO_1`, `KO_2`, and `KO_ctl`):

* `"results/bigwig/WT_1.bw"`
* `"results/bigwig/WT_2.bw"`
* `"results/bigwig/WT_ctl.bw"`
* `"results/bigwig/KO_1.bw"`
* `"results/bigwig/KO_2.bw"`
* `"results/bigwig/KO_ctl.bw"`

By explicitly naming files in this target rule, Snakemake is able to fill in the blanks for all wildcards encountered throughout the pipeline, making it so that the input and output file names match what is expected given this target rule.

