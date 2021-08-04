log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")
library("ggrepel")
library("tidyverse")
library("circlize")
library("EnhancedVolcano")
library("readr")
library("janitor")
library("GSA")
library("msigdbr")
library("fgsea")
library("shades")
library("biomaRt")

#Load contrast
elements_contrast <- snakemake@params[["contrast"]]
contrast <- paste(elements_contrast[1], "vs", elements_contrast[2], sep="-")

#Load gene set file names, species name, msigdb categories, descriptions option
elements <- snakemake@params[["gene_sets"]]
species <- snakemake@params[["species"]]
msigdb_categories <- snakemake@params[["msigdb_categories"]]
descriptions <- snakemake@params[["descriptions"]]

#Determine required columns 
if (descriptions == "yes") {
  required <- c("entrez_gene", "gs_name", "gs_description")
} else {
  required <- c("entrez_gene", "gs_name")
}

#Function to retrieve msigdbr species parameter 
if (any(msigdb_categories != "")) {
  if (species == "mouse") {
    species_msigdb <- "Mus musculus"
  } else if(species == "human") {
    species_msigdb <- "Homo sapiens"
  } else {
    stop("No msigdb entries for species: ",
         species)
  }
}


#Retrieve msigdb gene sets if specified
if (any(msigdb_categories != "")) {
  gene_sets <- lapply(msigdb_categories, msigdbr, species = species_msigdb) %>%
    bind_rows() %>%
    dplyr::select(all_of(required))
} else {
  gene_sets <- setNames(data.frame(matrix(ncol = length(required), nrow = 0)), required)
}

#Check to see if uploaded gene sets have required columns

##Function to check gene sets
check_gene_list <- function(filename, required_columns) {
  gene_set <-read.csv(filename)
  if (any(! required_columns %in% colnames(gene_set))) {
    missing <- required_columns[!required_columns %in% colnames(gene_set)] %>%
      paste(., collapse=",")
      
    stop(paste0("Dataframe of gene sets does not have ", missing, " column(s)"))
  }
  return(gene_set)
}

##If gene sets are uploaded, run check; if passes combine with earlier gene sets
if (any(elements != "")) {
  custom_gene_sets <- lapply(elements, check_gene_list, required) %>% 
    bind_rows() %>%
    dplyr::select(all_of(required))
  gene_sets <- gene_sets %>%
    rbind(custom_gene_sets)
}

gene_set_list <- split(x = gene_sets$entrez_gene, f = gene_sets$gs_name)

#Load full diffexp output and extract log2FoldChange
diff_exp <- read.csv(snakemake@input[['table']])%>% na.omit(cols=startsWith("padj"))
lfc <- dplyr::select(diff_exp, starts_with("gene_id"), starts_with("log2FoldChange"))

#Convert ensembl ID to entrez ID
species <- "mouse"
if (species == "mouse") {
	        ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
} else if(species == "human") {
	        ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
} else {
	        stop("No annotation (biomart) specified for organism:",
		                  species)
}

lfc_gene_ids  <- as.vector(lfc$gene_id)

#Get entrez ID's for our input genes and remove duplicates
lfc_names <- getBM(attributes=c('entrezgene_id', 'ensembl_gene_id'),
             values = lfc_gene_ids,
             mart = ensembl) %>%
  group_by(ensembl_gene_id) %>%
  filter(!any(row_number() > 1)) %>%
  na.omit()

#Select only the input genes that have corresponding entrez ID's
lfc_both_names <- lfc %>%
  subset(gene_id %in% lfc_names$ensembl_gene_id) %>%
  arrange(gene_id)

#Make sure all of the entrez ID's from biomaRt have ensembl ID's
lfc_names <- lfc_names %>%
  subset(ensembl_gene_id %in% lfc_both_names$gene_id) %>%
  arrange(ensembl_gene_id)

lfc_entrez_list <- cbind(lfc_both_names, lfc_names) %>%
  dplyr::select(starts_with("entrezgene_id"), starts_with("log2FoldChange")) %>%
  deframe()

#Run fgsea
fgseaRes <- fgsea(pathways = gene_set_list,
		  stats = lfc_entrez_list,
		  eps = 0.00,
		  minSize  = 15,
		  maxSize  = 500)

#Collect results
res <- fgseaRes[order(padj), ]

#If gene sets have descriptions, add to results
if (descriptions == "yes") {
  names_list <- as.vector(res$pathway)
  gene_des <- gene_sets %>%
	  dplyr::select(gs_name, gs_description) %>%
		distinct() %>%
    filter(gs_name %in% names_list) %>%
    dplyr::rename(pathway = gs_name)

	res <- list (res, gene_des) %>%
	  reduce(full_join, by = "pathway")
} 

#Save results as RDS object
saveRDS(res, file = snakemake@output[['rds']])

#Remove leadingEdge list column
res_no_le <- res %>%
	subset(select = -c(leadingEdge)) 

#Write results to csv
write_csv(res_no_le, snakemake@output[["table"]])

#Select the significant gene sets (if any) and add ensembl ID's for genes in those gene sets
if (any(res$padj < 0.05)) {
  unpack_le <- function(res_sig) {
    genes <- as.vector(res_sig$leadingEdge)
    return(genes)
  }
  get_BM_vector <- function(genes) {
    genes_ids_df <- getBM(attributes=c('entrezgene_id','ensembl_gene_id'),
                          filters = 'entrezgene_id',
                          values = genes,
                          mart = ensembl)
    genes_ids <- genes_ids_df$ensembl_gene_id
    return(genes_ids)
  }
  res_sig <- filter(res, padj < 0.05) 
  genes <- apply(res_sig, 1, unpack_le)
  genes_ids <- sapply(genes, get_BM_vector) 
  res_sig$leadingEdgeEnsembl <- genes_ids
  
  #Save significant results as RDS object
  saveRDS(res_sig, snakemake@output[["rds_sig"]]) 
  
  #Remove leadingEdge list column
  res_sig <- res_sig %>%
    subset(select = -c(leadingEdge, leadingEdgeEnsembl)) %>%
    arrange(padj)
  
  #Write significant results to csv
  write_csv(res_sig, snakemake@output[["table_sig"]])
}



#Visualize results
##Collect up to the top 5 pathways
if (nrow(res) >= 5) {
  n <- 5
} else {
  n <- nrow(res)
}
subset_pathway <- res[1:n,]
top_pathway_names <- subset_pathway$pathway

##Create enrichment plots
pdf(snakemake@output[["top_gene_sets_pdf"]])

enrichment_plots <- function(top_pathway_names) {
  plotEnrichment(pathway = gene_set_list[[top_pathway_names]], lfc_entrez_list) +
    labs(title = top_pathway_names) +
    theme(plot.title = element_text(hjust = 0.5, face="bold"))
} 

plt <- lapply(top_pathway_names,enrichment_plots)
print(plt)
dev.off()

#Create top pathway dotplot 
##Plot up to top 25 pathways and save to pdf
if (nrow(res) >= 25) {
  m <- 25
} else {
  m <- nrow(res)
}
selected_res <- res[1:m,] %>%
	  na.omit()

##Order results by NES
selected_res <- selected_res %>%
	arrange(NES) %>%
	mutate(pathway=fct_inorder(pathway))

gg <- ggplot(data=selected_res, aes(x=x, y=y)) + 
	geom_segment(aes(x=0, xend=NES, y=pathway, yend=pathway)) + 
	geom_point(aes(size=size, color= -log10(padj), x=NES, y=pathway)) + 
	theme_light() + 
	ggtitle("Enriched Pathways") + 
	lightness(scale_colour_distiller(palette = "YlGnBu", direction=1), scalefac(0.90)) + 
	xlab("Normalized Enrichment Score") + 
	ylab("Pathway") + 
	expand_limits(x=0) + 
	theme(
	        panel.grid.major.y = element_blank(),
		  panel.border = element_blank(),
		  axis.ticks.y = element_blank()
		  )

ggsave(plot=gg, scale = 2, filename = snakemake@output[["pathways_pdf"]], dpi = 600, width = 8, height = 6, units = "in")
