rule count_matrix:
    input:
        expand("star/{unit.sample}-{unit.unit}/ReadsPerGene.out.tab",
        unit=samples.itertuples())
    output:
        "counts/all.tsv"
    params:
        samples=samples["sample"].tolist(),
        strand=get_strandness(samples)
    resources:
        mem_mb = 10000
    conda:
       "../envs/pandas.yaml"
    log:
        "logs/counts/count_matrix.log"
    script:
        "../scripts/count-matrix.py"

rule deseq2_init:
    input:
        counts="counts/all.tsv"
    output:
        "deseq2/all.rds",
    params:
        samples=config["samples"],
        design=config["diffexp"]["design"],
        species=config["diffexp"]["organism"],
        collapse=config["diffexp"]["collapse"]
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/init.log"
    threads: get_deseq2_threads()
    script:
        "../scripts/deseq2-init.R"

rule pca:
    input:
        "deseq2/all.rds"
    output:
        pca_svg=report("results/pca.svg", caption = "../report/pca.rst"),
        pca_pdf="results/pca.pdf"
    params:
        fill=config["pca"]["fill"],
        color=config["pca"]["color"],
        labels=config["pca"]["labels"],
        batch=config["params"]["batch_adjust"]

    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/pca.log"
    script:
        "../scripts/plot-pca.R"

rule deseq2:
    input:
        "deseq2/all.rds"
    output:
        table= report("results/diffexp/{contrast}.diffexp.csv", caption = "../report/diffexp.rst"), 
        ma_plot=report("results/diffexp/{contrast}.ma-plot.svg", caption = "../report/ma.rst"),
        ma_pdf="results/diffexp/{contrast}.ma-plot.pdf",
        up=report("results/diffexp/deg-sig-up_{contrast}.csv", caption = "../report/up.rst"),
        down=report("results/diffexp/deg-sig-down_{contrast}.csv", caption = "../report/down.rst")
    params:
        contrast=get_contrast,
        annotationhub=config["diffexp"]["organism"]
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/{contrast}.diffexp.log"
    threads: get_deseq2_threads
    script:
        "../scripts/deseq2.R"

rule heatmap:
    input:
        "deseq2/all.rds",
        table ="results/diffexp/{contrast}.diffexp.csv",
        up = "results/diffexp/deg-sig-up_{contrast}.csv",
        down = "results/diffexp/deg-sig-down_{contrast}.csv"
    output:
        pdf=report("results/heatmap_{contrast}.pdf", caption = "../report/heatmap.rst"),
        counts_heatmap_data = "results/heatmap-app/data/hm_counts_sig_{contrast}.rds",
        lfc_heatmap_data = "results/heatmap-app/data/hm_fc_sig_{contrast}.rds",
        z_stat_heatmap_data = "results/heatmap-app/data/hm_z_sig_{contrast}.rds"
    params:
        contrast = get_contrast,
        species = config["diffexp"]["organism"],
        batch=config["params"]["batch_adjust"]
    threads: 1
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/heatmap_{contrast}.log"
    script:
        "../scripts/plot-heatmap.R"


if config["params"]["enrichment"]:
    rule fgsea:
        input:
            table ="results/diffexp/{contrast}.diffexp.csv"
        output:
            table= "results/fgsea/results_{contrast}.csv",
            top_gene_sets_pdf = report("results/fgsea/top_gene_sets_{contrast}.pdf", "../report/top_gene_sets.rst"),
            rds = "results/fgsea/res_{contrast}.rds",
            rds_sig = "results/fgsea/res_sig_{contrast}.rds",
            table_sig = "results/fgsea/results_sig_{contrast}.csv",
            pathways_pdf = report("results/fgsea/enriched_pathways_{contrast}.pdf", caption = "../report/top_pathways.rst"),
        params:
            contrast = get_contrast,
            gene_sets = config["fgsea"]["gene_sets"],
            species = config["diffexp"]["organism"],
            descriptions = config["fgsea"]["descriptions"],
            msigdb_categories = config["fgsea"]["msigdb_categories"]
        threads: 1
        resources:
            mem_mb = 60000
        conda:
            "../envs/deseq2.yaml"
        log:
            "logs/fgsea_{contrast}.log"
        script:
            "../scripts/fgsea.R"

rule volcano:
    input:
        table ="results/diffexp/{contrast}.diffexp.csv",
    output:
        pdf=report("results/volcano_{contrast}.pdf", caption = "../report/volcano.rst"),
    threads: 1
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/volcano_{contrast}.log"
    script:
        "../scripts/plot-volcano.R"
