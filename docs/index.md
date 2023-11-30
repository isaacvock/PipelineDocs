# Welcome my Snakemake pipeline deployment docs!

Biologists often rely on computational pipelines to process and analyze their data. A pipeline is simply a set of operations to perform on data so that you can extract useful insights from the data. With the omics-revolution came a spike in demand for these kinds of pipelines, as methods like RNA-seq generate massive raw datasets that require computationally intensive, multi-step processing to interpret. One of the best way to build such pipelines is with workflow managers, of which there are several popular and powerful options (Snakemake, NextFlow, Cromwell, etc.).

During my PhD, I became incredibly fond of Snakemake, and quickly realized just how transformative the modern generation of workflow managers are. They make building high quality, reproducible, and well organized pipelines easier than ever. For that reason, I devoted a lot of time to learning the ins and outs of Snakemake and the pipeline development process. In the process, I have developed a number of Snakemake pipelines for processing data relevant to my thesis project, and to the larger work of the [Simon lab](https://simonlab.yale.edu/). This website is a one-stop-shop for all documentation related to running these pipelines. In addition, I have attempted to generalize the instructions for running a Snakemake pipeline so that they apply equally well to all of the standardized Snakemake pipelines posted on the [Snakemake workflow catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=cbg-ethz/). Thus, I hope this a useful resource for anyone interested in using any of the pipelines I have developed as well as the treasure trove of amazing Snakemake pipelines developed by others.

## Where to go next

If you are here for information regarding one of my pipelines, I suggest checking out (in this order):

1. The relevant introduction page. The pipelines I developed whose documentation is currently hosted here include:
    - [bam2bakR](intro/bam2bakR_i.md)
    - [THE_Aligner](intro/aligner_i.md)
    - [PROseq_etal](intro/proseq_i.md)
1. The general pipeline deployment instructions [page](deploy.md), or the Yale HPC-specific one if [relevant](simon.md). The latter also provides some information that may be useful to anyone running these pipelines on a shared computing system in which it is possible to request multiple jobs in parallel.
1. The relevant configuration page:
    - [bam2bakR](configuration/bam2bakR_c.md)
    - [THE_Aligner](configuration/aligner_c.md)
    - [PROseq_etal](configuration/proseq_c.md)
1. The tips and tricks [page](pragmatism.md) for general advice about Snakemake which might help demystify these pipelines a bit.
1. Information about the output produced by each pipeline:
    - [bam2bakR](output/bam2bakR_o.md)
    - [THE_Aligner](output/aligner_o.md)
    - [PROseq_etal](output/proseq_o.md)

While not currently populated (for the most part), other pages that may be useful to check out in the future if you are using my pipelines are:

1. The FAQs for each pipeline, with some pipeline-specific debugging tips.
    - [bam2bakR's](faqs/bam2bakR_f.md) already has some useful information there.
1. The architecture section, which will go into great detail about how each and every step of the pipeline works.


Everyone else should head to the general pipeline deployment instructions [page](deploy.md), or the Yale HPC-specific one if [relevant](simon.md), and then to the tips and tricks [page](pragmatism.md)!