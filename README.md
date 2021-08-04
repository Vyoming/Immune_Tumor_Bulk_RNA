# RNAseq Snakemake workflow

This workflow performs a differential expression analysis using [STAR](https://github.com/alexdobin/STAR) and [DESeq2](http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html). It performs a wide range of quality control steps and
bundles there results in a qc report via [MultiQC](https://multiqc.info/). Reported results include:
* PCA plot of all samples
* differentially up and down regulated genes for each analysis (ie treatment vs control)
* MA plot for each analysis
* heatmap of counts and log2 Fold Changes
* lollipop plot of enriched gene sets

## Usage

#### Step 1: Configure workflow

You will likely run these analyses on the HPCC (currently Elzar), but you will want access to the results locally.
Our strategy is to have the analysis repository on the local machine and on the HPCC, in the same file path location
(relative to the home directory). 

We will start by creating the local repository (on your computer): 
* Click on 'Use this template', which will copy the content of this repository to a new remote repositiory on github.
Find a good name for your analysis and use it as the name for the repository;
* Copy git link from this new repository, move to the parent folder for your RNA-seq analysis (ie where you want to save your new analyses/results) and
create a local repository on your computer using git clone;
* Move into the newly created directory (repository);
* Configure the workflow by editing the file `config.yaml`;
    + Choose mouse or human for the STAR Index;
    + Modify the pca columns, design for deseq2, organism of interest, the contrasts, gene sets for (optional) enrichment analysis;
    + Add gene set files to gene_sets directory
* Configure the `samples.txt` and `units.txt` files with project specific annotations and file;
    + Hint: Check for rogue leading/trailing spaces in your file paths/sample names etc;
* Git push changes to master.

We will now set up the analysis on the HPCC:
* Copy the git link again and create a repository on the cluster (same path, see above) using git clone;
* Move into the newly created directory (repository);
* Create a reads directory and copy over the reads (fastq.gz files) from the sequencing facility.

To run the analysis, activate your conda snakemake environment:
```
    conda activate snakemake
```

#### Step 2: Execute workflow

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

Execute the workflow locally via

    snakemake --use-conda --cores $N

using `$N` cores or run it on the [UGE cluster environment](https://github.com/meyer-lab-cshl/snakemake-uge)

    snakemake --use-conda --profile uge

#### Step 3: Investigate results

After successful execution, you can create a self-contained interactive HTML report via:

    snakemake --report report.html

* to look at the results, use secure copy `scp` to transfer the `report.html` and the quality control report `qc/multiqc.html` from the HPCC to your
local computer (this is where storing them in the same location relative to the home directory comes in handy);
* explore the results; 
* forward to your collaborators.

This workflow is based on [rna-seq-star-deseq2](https://github.com/snakemake-workflows/rna-seq-star-deseq2/releases).
