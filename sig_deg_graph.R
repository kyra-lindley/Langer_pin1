# This file takes the results of the the deseq_init and makes a heatmap looking at the top DEG

# To load the data 
res <- readRDS("/Users/lindleyk/Langer/Bulk/deseq/res_PSC_Pin1_TGFBvsUN_rld.rds")

topGenes <- head(order(res$padj),50)

df <- as.data.frame(colData(rld))
if (length(subset_cols)==1) {
  annot <- as.data.frame(cbind(rownames(df), paste(df[[subset_cols[1]]])))
  names(annot) <- c("SampleID", subset_cols[1])
  rownames(annot) <- annot$SampleID
  annot$SampleID <- NULL
} else {
  annot <- df[,subset_cols]
}

sig_genes <- res %>%
  as.data.frame() %>%
  rownames_to_column("ensembl") %>%
  filter(!is.na(padj), padj < 0.05) %>%
  arrange(padj)

topGenes <- sig_genes$ensembl
# Replace NAs with ENSEMBL IDs
gene_symbols[is.na(gene_symbols)] <- topGenes[is.na(gene_symbols)]
df <- as.data.frame(colData(rld))
annot <- df[, c("group")]  # Or whatever metadata you want
# Convert to data frame and assign sample names as rownames
annot <- data.frame(Group = rld$group)
rownames(annot) <- colnames(rld)
# If your topGenes have versions (e.g. ENSG00000123456.1), strip them:
topGenes_clean <- gsub("\\..*", "", topGenes)

# Now subset gene_symbols to just those genes, preserving order:
gene_labels <- gene_symbols[topGenes_clean]

# Fill in missing symbols with ENSEMBL IDs so there are no blanks:
gene_labels[is.na(gene_labels)] <- topGenes_clean[is.na(gene_labels)]

pheatmap(
  assay(rld)[topGenes, ],
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "row",
  show_rownames = TRUE,
  labels_row = gene_labels,
  annotation_col = annot,
  main = paste("Heatmap of statistically significant DE genes")
)
