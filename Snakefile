import pandas as pd
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.14.0")


##### load config and sample sheets #####
configfile: "config.yaml"
validate(config, schema="schemas/config.schema.yaml")

samples = pd.read_table(config['samples']).set_index(["sample", "unit"],
    drop=False)
validate(samples, schema="schemas/samples.schema.yaml")

# enforce str in index
samples.index = samples.index.set_levels([i.astype(str) for i in samples.index.levels])

##### target rules #####
if config["params"]["enrichment"] == "yes":
  rule all:
        input:
            expand(["results/diffexp/{contrast}.diffexp.csv",
                    "results/heatmap_{contrast}.pdf", 
                    "results/fgsea/results_{contrast}.csv",
                    "results/volcano_{contrast}.pdf"],
                    contrast=config["diffexp"]["contrasts"]),
            "results/pca.svg",
            "qc/multiqc_report.html",

else:
    rule all:
        input:
            expand("star/{samples.sample}-{samples.unit}/Aligned.sortedByCoord.out.bam.bai",
                samples=samples.itertuples()),
            expand(["results/diffexp/{contrast}.diffexp.csv",
                    "results/diffexp/{contrast}.ma-plot.svg",
                    "results/heatmap_{contrast}.pdf", 
                    "results/heatmap-app/data/hm_counts_sig_{contrast}.rds",
                    "results/heatmap-app/data/hm_fc_sig_{contrast}.rds",
                    "results/heatmap-app/data/hm_z_sig_{contrast}.rds",
                    "results/volcano_{contrast}.pdf"],
                    contrast=config["diffexp"]["contrasts"]),
            "results/pca.svg",
            "qc/multiqc_report.html",

##### setup report #####
report: "report/workflow.rst"

##### load rules #####
include: "rules/functions.smk"
include: "rules/trim.smk"
include: "rules/align.smk"
include: "rules/diffexp.smk"
include: "rules/qc.smk"
