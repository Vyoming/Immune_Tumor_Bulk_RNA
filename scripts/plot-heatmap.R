log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

## libraries ####
library("DESeq2")
library("cowplot")
library("limma")
library("ggrepel")
library("tidyverse")
library("scales")
library("ggnewscale")
library("ComplexHeatmap")
library("circlize")
library("EnhancedVolcano")
library("readr")
library("janitor")
library("biomaRt")

## functions ####
theme_set(theme_half_open(font_family = 'sans' ))

title_size <- 10
text_size <- 8
point_size <- 0.8

## settings ####
map_color_range <- function(matrix, color_vec) {
    matrix %>%
        range(., na.rm=TRUE) %>%
        abs %>%
        max %>%
        seq(-., ., length.out = length(color_vec)) %>%
        colorRamp2(., color_vec)
}

## data ####
dds <- readRDS(snakemake@input[[1]])
elements <- snakemake@params[["contrast"]]
contrast <- paste(elements[1], "vs", elements[2], sep="-")
adjust_batch <- snakemake@params[['batch']]


## analysis ####
# keep contrast specific samples
dds <- dds[,colData(dds)$condition %in% elements]

# obtain normalized counts
counts <- rlog(dds, blind=FALSE)

if (adjust_batch != "") {
    if (!adjust_batch %in% colnames(colData(counts))) {
        stop("Column selected for batch correction prior to PCA not present")
    }
    assay(counts) <- limma::removeBatchEffect(assay(counts),
                                              colData(counts)[, adjust_batch])
}

# Reformat and scale counts
counts_scaled <- t(apply(assay(counts), 1, scale))
colnames(counts_scaled) <- colnames(assay(counts))

# Load up and downregulated genes and extract gene list
up_reg <- read.csv(snakemake@input[['up']])
down_reg <- read.csv(snakemake@input[['down']])
sig_genes <- rbind(up_reg, down_reg)

fc <- sig_genes %>%
  dplyr::select(gene_id, gene_name, log2FoldChange)

se <- sig_genes %>%
  dplyr::select(gene_id, lfcSE)

# Extract counts and fc values for significant genes and order them by gene_id
combined <- counts_scaled %>%
  as.data.frame %>%
  rownames_to_column(var = "gene_id") %>%
  inner_join(fc, by="gene_id") %>%
  inner_join(se, by="gene_id") %>%
  dplyr::select(gene_id, gene_name, log2FoldChange, lfcSE,everything()) %>%
  arrange(gene_id) %>%
  as.data.frame

# Set heatmap categories and condition type
if (!"labels" %in% colnames(colData(dds))) {
    colData(dds)$labels <- colData(dds)$condition
}

# Set heatmap colors, themes, and sizes
cond_type_vals <- c('#d95f02','#7570b3')
names(cond_type_vals) <- unique(colData(dds)$labels)
color <- list(cond_type = cond_type_vals)
counts_color <- c('#a6611a','#dfc27d','#f5f5f5','#80cdc1','#018571')
fc_color <- c('#0571b0','#92c5de','#f7f7f7','#f4a582','#ca0020')
z_stat_color <- c('#e66101','#fdb863','#f7f7f7','#b2abd2','#5e3c99')

# Create column annotations for rnaseq data
coldata <- colData(dds) %>%
  as_tibble

columnAnno <- columnAnnotation(
                   cond_type = coldata$labels,
                   col = list(cond_type = color$cond_type),
                   annotation_legend_param = list(cond_type = list(title = "Condition")),
                   show_annotation_name = FALSE,
                   annotation_name_gp = gpar(fontsize = text_size)
                   )

# Create the heatmap
counts_matrix <- combined %>%
  dplyr::select(-c(gene_id,  log2FoldChange, lfcSE)) %>%
  column_to_rownames("gene_name") %>%
  as.matrix

fc_matrix <- combined %>%
  dplyr::select(gene_name, log2FoldChange) %>%
  column_to_rownames("gene_name") %>%
  as.matrix

z_stat_matrix <- combined %>%
  mutate(z_stat = log2FoldChange / lfcSE) %>%
  dplyr::select(gene_name, z_stat) %>%
  column_to_rownames("gene_name") %>%
  as.matrix

hm_counts_sig <- Heatmap(counts_matrix,
                         name="Counts",
                         col=map_color_range(counts_matrix, counts_color),
                         cluster_rows=TRUE,
                         cluster_columns=FALSE,
                         show_row_dend = FALSE,
                         show_row_names=FALSE,
                         show_column_names=TRUE,
                         top_annotation = columnAnno,
                         row_names_gp = gpar(fontsize = text_size),
                         column_names_gp = gpar(fontsize = text_size),
                         heatmap_legend_param =
                                 list(direction = "horizontal",
                                 title_gp = gpar(fontsize = text_size,
                                                      fontface = "bold"),
                                      labels_gp = gpar(fontsize = text_size)),
                         use_raster=TRUE,
                         row_title=NULL,
                         column_title=NULL)

hm_fc <- Heatmap(fc_matrix,
         col=map_color_range(fc_matrix, fc_color),
         use_raster=TRUE,
         show_row_names=FALSE,
         show_row_dend = FALSE,
         show_column_names = FALSE,
                 row_title = NULL,
         column_title = NULL,
         cluster_columns=FALSE,
         cluster_rows=TRUE,
         row_names_gp = gpar(fontsize = text_size),
         column_names_gp = gpar(fontsize = text_size),
         heatmap_legend_param =
             list(direction = "horizontal",
                  title_gp = gpar(fontsize = text_size,
                          fontface = "bold"),
                  labels_gp = gpar(fontsize = text_size)),
         name = "log2(FC)",
         width = unit(3, "cm"))

hm_z <- Heatmap(z_stat_matrix,
                 col=map_color_range(z_stat_matrix, z_stat_color),
                 use_raster=TRUE,
                 show_row_names=FALSE,
                 show_row_dend = FALSE,
                 show_column_names = FALSE,
                 row_title = NULL,
                 column_title = NULL,
                 cluster_columns=FALSE,
                 cluster_rows=TRUE,
                 row_names_gp = gpar(fontsize = text_size),
                 column_names_gp = gpar(fontsize = text_size),
                 heatmap_legend_param =
                   list(direction = "horizontal",
                        title_gp = gpar(fontsize = text_size,
                                        fontface = "bold"),
                        labels_gp = gpar(fontsize = text_size)),
                 name = "z_statistic",
                 width = unit(3, "cm"))

pdf(snakemake@output[['pdf']], width=190/25.4, height =150/25.4)
draw(hm_counts_sig + hm_fc, heatmap_legend_side = "bottom",
          merge_legend = TRUE)
dev.off()

#Save counts heatmap
saveRDS(hm_counts_sig, file = snakemake@output[['counts_heatmap_data']])
saveRDS(hm_fc, file = snakemake@output[['lfc_heatmap_data']])
saveRDS(hm_z, file = snakemake@output[['z-stat_heatmap_data']])
