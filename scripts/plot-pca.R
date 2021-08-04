log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")
library("ggplot2")
library("cowplot")
library("limma")
library("ggrepel")

dds <- readRDS(snakemake@input[[1]])
pca_color <- snakemake@params[['color']]
pca_fill <- snakemake@params[['fill']]
pca_labels <- snakemake@params[['labels']]
pca_adjust_batch <- snakemake@params[['batch']]


if (all(c(pca_color, pca_fill) != "")) {
    intgroup <- c(pca_color, pca_fill)
} else if (pca_color != "") {
    intgroup <- pca_color
} else if (pca_fill != "") {
    intgroup <- pca_fill
} else {
    stop("At least one of fill or color have to be specified")
}

if (pca_labels != "") {
    if (!pca_labels %in% intgroup) {
        intgroup <- c(intgroup, pca_labels)
    }
}

head(colnames(colData(dds)))
if (!all(intgroup %in% colnames(colData(dds)))) {
    missing <- intgroup[!intgroup %in% colnames(colData(dds))]
    stop(paste(missing, collapse=","),
         "columns selected for PCA highlighting are missing")
}

# obtain normalized counts
counts <- rlog(dds, blind=FALSE)

if (pca_adjust_batch != "") {
    if (! pca_adjust_batch %in% colnames(colData(counts))) {
        stop("Column selected for batch correction prior to PCA not present")
    }
    assay(counts) <- limma::removeBatchEffect(assay(counts),
                                              colData(counts)[, pca_adjust_batch])

}
pcaData <- plotPCA(counts, intgroup = intgroup, returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))


color_sym <- sym(pca_color)
p <- ggplot(pcaData, aes(x = PC1, y = PC2))

if (all(c(pca_color, pca_fill) != "")) {
    fill_sym <- sym(pca_fill)
    p <- p +
        geom_point(aes(color = !!color_sym, fill = !!fill_sym),
                   pch=22,  size=3, stroke=2)
} else {
    p <- p + geom_point(aes(color = !!color_sym))
}

if (pca_labels != "") {
    label_sym <- sym(pca_labels)
    p <- p + geom_text_repel(aes(label=!!label_sym), size = 3, color="black",
                             segment.size = 0.3, segment.color = "black",
                             min.segment.length = 0.2, point.padding = 0.3)
}

p <- p +
    scale_color_brewer(type="qual", palette="Dark2") +
    labs(x=paste0("PC1: ", percentVar[1], "% variance"),
         y=paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed() +
    theme_cowplot() +
    theme(legend.position = "bottom",
          legend.justification = 0)
ggsave(plot=p, height=4.5, width=7.5, filename = "results/pca.pdf")
ggsave(plot=p, height=4.5, width=7.5, filename = "results/pca.svg")


