## bam2bakR FAQs

### Why are there a low number of T-to-C mutations in the cB?

If there are very few T-to-C mutations in the final cB.csv file (e.g., if sample-wide mutation rates in +s4U samples are < 0.003), then you may have used the incorrect value for the `strandedness` parameter in the config. One way to tell if this is the case is by looking at one of the +s4U sample counts.csv files in `results/counts/` and checking for an abundance of A-to-G mutations. If this is the case, flip the value of `strandedness` to the opposite of whatever you used.

Related to the first point, a good sanity check after running the pipeline is going into R and checking the raw mutation rates as such:

```r
library(data.table)

# To unzip and read cB, also need to have R.utils package installed
cB <- fread("path/to/cB.csv.gz")

# Assess sample-wide T-to-C mutation rate in each sample
cB[,.(mutrate = sum(TC*n)/sum(nT*n), by = sample]
  # Want to see that +s4U samples has higher mutation rate than -s4U samples
```

Similarly, checking a counts.csv file for an abundance of A-to-G mutations can be done as follows:

```r
library(data.table)

counts <- fread("path/to/+s4U/counts.csv.gz")

## Check if A-to-G mutation rate is higher than T-to-C mutation rate:

# A-to-G mutation rate
sum(counts$AG)/sum(counts$nA)

# T-to-C mutation rate
sum(counts$TC)/sum(counts$nT)
```

### How can I run older versions of bam2bakR?

There are a number of previous versions of bam2bakR that have "tags" on Github. If you go to the Github and click on the box that says "main" (which represents the branch being displayed), a dropdown menu will appear. In that menu, you will see two subsections labeled Branches and Tags. Click Tags to see the available tags. 

Using one of the tagged versions is easy if you deployed bam2bakR using Snakedeploy. To do so, navigate to the Snakefile located in the "workflow" director created when you deploy the pipeline. It will look something like:

```
from snakemake.utils import min_version


min_version("6.10.0")


configfile: "config/config.yaml"


# declare https://github.com/simonlabcode/bam2bakR as a module
module bam2bakR:
    snakefile:
        github("simonlabcode/bam2bakR", path="workflow/Snakefile", branch="main")
    config:
        config


# use all rules from https://github.com/simonlabcode/bam2bakR

```

The key part is:

```
# declare https://github.com/simonlabcode/bam2bakR as a module
module bam2bakR:
    snakefile:
        github("simonlabcode/bam2bakR", path="workflow/Snakefile", branch="main")
    config:
        config
```

To use a tagged version, change `branch="main"` to `tag="<name of tag>"`. For example, to revert to version 2.0.1:

```
# declare https://github.com/simonlabcode/bam2bakR as a module
module bam2bakR:
    snakefile:
        github("simonlabcode/bam2bakR", path="workflow/Snakefile", tag="2.0.1")
    config:
        config
```

This can also be set when deploying the workflow. For example, if you change your call to Snakedeploy to:

```
snakedeploy deploy-workflow https://github.com/simonlabcode/bam2bakR . --tag v2.0.1
```

this will deploy version 2.0.1.

One challenge you may run into when trying to revert to an older version is that the config file parameters in the newer version may differ from those used by the older version. In this case, when you run the pipeline, you will get `KeyErrors` when running Snakemake. You can use these to track down which config file parameters are needed. What is probably easier though is to go to the bam2bakR Github and get the config file from the relevant tag.


### Lots of NA XF or GF

If a large number of reads in the cB.csv.gz file have an XF/GF (i.e., the annotated feature the read was assigned to) value of NA (or __ambiguous in version 2.0.1 or earlier), this likely has a similar cause to the low T-to-C mutation problem discussed above. That is, the strandedness of the library may be improperly defined in the config. Thus, check that the `strandedness` parameter in the config file is properly set, and confirm the strandedness of your library.

Version 1.0.1 of bam2bakR and earlier had a bug that caused quantification of features to always assume forward strandedness. This bug was corrected in version 1.0.2. Thus, one solution is to rerun bam2bakR with version 1.0.2+, with the proper library strandedness noted in the config. If you have reverse stranded data that you would like to run through an older version of bam2bakR, a simple hack is to just flip the fastq files during alignment. That is, provide read 2 as the read 1 fastq and vice versa. This will reverse the strandedness of your library and make it effectively forward stranded.

