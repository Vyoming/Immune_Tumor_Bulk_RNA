library("ggrepel")
library("tidyverse")

#Load full diffexp output and extract log2FoldChange and padj
diffexp <- read.csv(snakemake@input[['table']])
fc_p_adj <- dplyr::select(diffexp, starts_with("gene_id"), starts_with("log2FoldChange"), starts_with("padj"), starts_with("gene_name")) %>%
  na.omit %>%
  rename(
	 log2FoldChange = starts_with("log2FoldChange"),
	 padj = starts_with("padj")
	 )

# The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
# Add a column to the data frame to specify if they are UP- or DOWN- regulated

# add a column of NAs
fc_p_adj$diff_expressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
fc_p_adj$diff_expressed[fc_p_adj$log2FoldChange > 0.6 & fc_p_adj$pad < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
fc_p_adj$diff_expressed[fc_p_adj$log2FoldChange < -0.6 & fc_p_adj$padj < 0.05] <- "DOWN"

# Now write down the name of genes beside the points...
# Create a new column "delabel" to the dataframe, that will contain the name of genes differentially expressed (NA in case they are not)
fc_p_adj$delabel <- NA
fc_p_adj$delabel[fc_p_adj$diff_expressed != "NO"] <- fc_p_adj$gene_name[fc_p_adj$diff_expressed != "NO"]

#Create and save plot
p <- ggplot(data=fc_p_adj, aes(x=log2FoldChange, y=-log10(padj), col=diff_expressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  ggtitle("-log10(p) vs log2FoldChange") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(plot=p, scale = 1, width = 15.4, height = 8.6, filename = snakemake@output[['pdf']])
