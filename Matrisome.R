# This file makes the comparison for the naba matrisome list, comparing DEGs and projecting the overlap, getting the z-score, and understanidng which genes are up or down regulated. 

res <- readRDS("/Users/lindleyk/Langer/Bulk/deseq/psc_scr_tgfb_vs_un/res_psc_scr_tgfb_vs_un.rds")
rld <- readRDS("/Users/lindleyk/Langer/Bulk/deseq/psc_scr_tgfb_vs_un/rlog_psc_scr_tgfb_vs_un.rds")
library(org.Hs.eg.db)
library(dplyr)
library(tibble)

res_df <- as.data.frame(res) %>%
  rownames_to_column("ensembl") %>%
  mutate(ensembl_clean = gsub("\\..*", "", ensembl)) %>%
  filter(!is.na(padj), padj < 0.05) %>%
  mutate(symbol = mapIds(org.Hs.eg.db,
                         keys = ensembl_clean,
                         column = "SYMBOL",
                         keytype = "ENSEMBL",
                         multiVals = "first")) %>%
  filter(!is.na(symbol))

library(GSEABase)

# Replace with your actual path to the GMT file
gmt_file <- "/Users/lindleyk/Desktop/NABA_MATRISOME.v2024.1.Hs.gmt"  # Example

gene_sets <- getGmt(gmt_file)

# Convert to a list for easier access
gene_sets_list <- geneIds(gene_sets)

# View available gene set names
names(gene_sets_list)[grepl("NABA", names(gene_sets_list))]

# Extract NABA_MATRISOME specifically
naba_matrisome <- gene_sets_list[["NABA_MATRISOME"]]

matrisome_DEGs <- res_df %>%
  filter(symbol %in% naba_matrisome)

rld_mat <- assay(rld)

# Subset rows matching Ensembl IDs (cleaned) from DEGs
heat_mat <- rld_mat[rownames(rld_mat) %in% matrisome_DEGs$ensembl_clean, ]

# Replace Ensembl rownames with gene symbols
rownames(heat_mat) <- matrisome_DEGs$symbol[match(rownames(heat_mat), matrisome_DEGs$ensembl_clean)]

# If you have a sample mapping vector and specific column order:
colnames(heat_mat)[colnames(heat_mat) %in% names(sample_map)] <- sample_map[colnames(heat_mat)[colnames(heat_mat) %in% names(sample_map)]]
heat_mat <- heat_mat[, desired_order]

library(pheatmap)

# === Step 1: Prepare heatmap matrix ===
rld_mat <- assay(rld)

heat_mat <- rld_mat[rownames(rld_mat) %in% matrisome_DEGs$ensembl_clean, ]
rownames(heat_mat) <- matrisome_DEGs$symbol[match(rownames(heat_mat), matrisome_DEGs$ensembl_clean)]

# Optional: Apply sample renaming and desired order
colnames(heat_mat)[colnames(heat_mat) %in% names(sample_map)] <- sample_map[colnames(heat_mat)[colnames(heat_mat) %in% names(sample_map)]]
heat_mat <- heat_mat[, desired_order]

# === Step 2: Create cleaned annotation_col ===
# Extract group names (e.g., shSCR_UN) by removing sample numbers
group_labels <- gsub(" [0-9]+$", "", colnames(heat_mat))        # e.g. "shSCR UN 1" → "shSCR UN"
group_labels <- gsub(" ", "_", group_labels)                    # e.g. "shSCR UN" → "shSCR_UN"
group_labels <- trimws(group_labels)

annotation_col <- data.frame(group = factor(group_labels,
                                            levels = c("shSCR_UN", "shSCR_TGFB", "shPIN1_UN", "shPIN1_TGFB")))
rownames(annotation_col) <- colnames(heat_mat)

# === Step 3: Define colors ===
ann_colors <- list(
  group = c(
    "shSCR_UN" = "black",
    "shSCR_TGFB" = "gray",
    "shPIN1_UN" = "red",
    "shPIN1_TGFB" = "pink"
  )
)

# === Step 4: Plot heatmap ===
pheatmap(
  heat_mat,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  fontsize_col = 8,
  main = "NABA Matrisome Genes Differentially Expressed"
)


library(dplyr)
library(ggplot2)
library(writexl)

# === Step 1: Join expression matrix with log2FC ===
gene_info <- matrisome_DEGs[, c("symbol", "log2FoldChange")]
heat_df <- as.data.frame(heat_mat)
heat_df$symbol <- rownames(heat_mat)
heat_df_joined <- inner_join(heat_df, gene_info, by = "symbol")

# === Step 2: Split into upregulated and downregulated subsets ===
heat_up <- heat_df_joined %>% filter(log2FoldChange > 0) %>% select(-log2FoldChange)
heat_down <- heat_df_joined %>% filter(log2FoldChange < 0) %>% select(-log2FoldChange)

# === Step 3: Z-score across genes ===
zscore_matrix <- function(mat) {
  t(scale(t(as.matrix(mat[, -ncol(mat)]))))  # exclude symbol column
}

z_up <- zscore_matrix(heat_up)
z_down <- zscore_matrix(heat_down)

# === Step 4: Average z-score per sample ===
zscore_summary <- data.frame(
  Sample = rep(colnames(z_up), 2),
  Zscore = c(colMeans(z_up, na.rm = TRUE), colMeans(z_down, na.rm = TRUE)),
  Direction = rep(c("NABA-up", "NABA-down"), each = ncol(z_up))
)

# === Step 5: Set sample order for plotting ===
zscore_summary$Sample <- factor(zscore_summary$Sample, levels = desired_order)

# === Step 6: Bar plot of cumulative z-scores ===
ggplot(zscore_summary, aes(x = Sample, y = Zscore, fill = Direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Cumulative Z-Score of NABA Matrisome-Regulated Genes",
       y = "Average Z-Score", x = NULL) +
  scale_fill_manual(values = c("NABA-up" = "firebrick", "NABA-down" = "steelblue")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# === Step 7: Export z-score summary (optional) ===
write_xlsx(zscore_summary, path = "/Users/lindleyk/Desktop/shSCR_TGFB_vs_UN_zscore_summary_NABA_matrisome.xlsx")

# === Step 8: Count NABA matrisome DEGs per sample ===
count_up <- colSums(!is.na(heat_up[, -ncol(heat_up)]))    # exclude symbol
count_down <- colSums(!is.na(heat_down[, -ncol(heat_down)]))

# === Step 9: Create summary table ===
matrisome_direction_summary <- tibble(
  Sample = rep(colnames(heat_mat), 2),
  Direction = rep(c("NABA-up", "NABA-down"), each = ncol(heat_mat)),
  GeneCount = c(count_up, count_down)
)

# === Step 10: Total unique DEGs in each direction ===
N_up <- nrow(heat_up)
N_down <- nrow(heat_down)

# === Step 11: Print summary tables ===
print(matrisome_direction_summary)

print(tibble(
  Direction = c("NABA-up", "NABA-down"),
  Total_NABA_DEGs = c(N_up, N_down)
))
