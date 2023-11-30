## Alternative Strategies for Running Snakemake Pipelines

An alternative way to run a Snakemake pipeline is to clone the pipeline's repo locally. This makes updating the pipeline a bit of a hassle, as it is non-trivial to pull any updates from the pipeline's repo if you edited the config file. That being said, cloning the full repo can be useful if you need to make changes to the workflow due to idiosyncracies in your particular data. Therefore, this section discusses how to run Snakemake pipelines in this way.

First, clone the pipeline's repository to wherever you would like on your system. You will eventually be navigating to this repo directory in the terminal and running Snakemake from inside the directory, so make sure your chosen location is conducive to this. Navigate to the directory in the terminal and run:

``` bash
$ git clone <path to pipeline .git>
$ cd <pipeline name>
```
You should be in the pipeline's repo directory now!

Next, make sure you have mamba/conda installed, as described in Step 2 of the [Snakedeploy route](deploy.md) (i.e., installing mamba/conda).

Pipeline's I developed as well as standardized pipelines posted on the [Snakemake workflow catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=cbg-ethz/) are compatible with Snakemake's [--use-conda option](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html). This will cause Snakemake to automatically create and activate conda environments for each step of the workflow to run inside. If you want to use this functionality, you can start by creating a simple conda environment that contains snakemake, as such:

``` bash
mamba create -c conda-forge -c bioconda --name snakemake snakemake
```

You would then run the pipeline with the `snakemake` environment activated with:

``` bash
snakemake cores all --use-conda
```

where `cores all` is a convenient way to tell Snakemake to make use of all available cpus (all can be replaced with an explicit number as was shown in the installation/pipeline running instructions above).