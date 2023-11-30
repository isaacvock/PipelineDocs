## Alternative Strategies for Running THE_Aligner

An alternative way to run THE_Aligner is to clone the THE_Aligner repo locally. This makes updating THE_Aligner a bit of a hassle, as it is non-trivial to pull any updates from the THE_Aligner repo if you edited the config file. That being said, cloning the full repo can be useful if you need to make changes to the workflow due to idiosyncracies in your particular data. Therefore, this section discusses how to run THE_Aligner in this way.

Clone the THE_Aligner repository to wherever you would like on your system. You will eventually be navigating to this repo directory in the terminal and running Snakemake from inside the directory, so make sure your chosen location is conducive to this. Navigate to the directory in the terminal and run:

``` bash
$ git clone https://github.com/isaacvock/THE_Aligner.git
$ cd THE_Aligner
```
You should be in the THE_Aligner repo directory now!

Step 1 of the Snakedeploy route (i.e., installing mamba/conda) is also the same for the Activating --use-conda and Running all rules inside one environment routes.

THE_Aligner is compatible with Snakemake's [--use-conda option](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html). This will cause Snakemake to automatically create and activate conda environments for each step of the workflow to run inside. If you want to use this functionality, you can start by creating a simple conda environment that contains snakemake, as such:

``` bash
mamba create -c conda-forge -c bioconda --name snakemake snakemake
```

You would then run the pipeline with the `snakemake` environment activated with:

``` bash
snakemake cores all --use-conda
```

where `cores all` is a convenient way to tell Snakemake to make use of all available cpus (all can be replaced with an explicit number as was shown in the installation/pipeline running instructions above).