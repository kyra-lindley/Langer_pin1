---
title: "deseq_init"
output: html_document
date: "2025-05-19"
---
This file was used to make deseq objects for further analysis. The input was metadata and counts file. Output was des, res, rld.

knitr::opts_chunk$set(echo = TRUE)
library(dplyr)
library(DESeq2)
library(readxl)
library(AnnotationDbi)
library(org.Hs.eg.db)

# === USER INPUTS === #
counts_file <- "/Users/lindleyk/Langer/Bulk/counts.txt"
metadata_file <- "/Users/lindleyk/Langer/Bulk/metadata_deseq2.xlsx"
sampleID <- "Sample"
group_col <- "group"
contrast = c("group", "CAFPIN_UN", "CAFSCR_UN")
# c(target, baseline)
output_rds <-"/Users/lindleyk/Langer/Bulk/deseq/CAF/shPin1_un_vs_shSCR_un/dds_CAF_shPin1_un_vs_shSCR_un.rds" 
output_rld <- "/Users/lindleyk/Langer/Bulk/deseq/CAF/shPin1_un_vs_shSCR_un/rlog_shPin1_un_vs_shSCR_un_rld.rds"
output_res <- "/Users/lindleyk/Langer/Bulk/deseq/CAF/shPin1_un_vs_shSCR_un/res_shPin1_un_vs_shSCR_un_rld.rds"

# === LOAD DATA === #
md <- read_excel(metadata_file)
md[[sampleID]] <- gsub("_Aligned.*", "", basename(md[[sampleID]]))
md <- md[order(md[[sampleID]]), ]

counts <- read.table(counts_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
colnames(counts) <- gsub("_Aligned.*", "", basename(colnames(counts)))
counts <- counts[, order(colnames(counts))]

# Keep only samples you're interested in
#keep_cols <- c("PSCPINT1", "PSCPINT2", "PSCPINT3", "PSCPINU1", "PSCPINU2", "PSCPINU3","PSCSCRT1", "PSCSCRT2", "PSCSCRT3", "PSCSCRU1", "PSCSCRU2", "PSCSCRU3")
keep_cols <-c("CAFPINT1", "CAFPINT2", "CAFPINT3", "CAFPINU1", "CAFPINU2", "CAFPINU3", "CAFSCRT1", "CAFSCRT2", "CAFSCRT3", "CAFSCRU1","CAFSCRU2", "CAFSCRU3")
counts <- counts[, colnames(counts) %in% keep_cols]

md <- md[md[[sampleID]] %in% keep_cols, ]
rownames(md) <- md[[sampleID]]

# Double-check alignment
stopifnot(all(rownames(md) == colnames(counts)))

# === DESeq2 SETUP === #
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = md,
                              design = as.formula(paste("~", group_col)))
dds <- dds[rowSums(counts(dds)) >= 1, ]
dds <- DESeq(dds)

# Ensembl -> Gene symbols
ensembl_ids <- gsub("\\..*", "", rownames(dds))
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = ensembl_ids,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)
rowData(dds)$symbol <- gene_symbols
head(rowData(dds)$symbol)

# Save DESeq object
saveRDS(dds, file = output_rds)

# === RLOG NORMALIZATION === #
rld <- rlog(dds, blind = FALSE)
saveRDS(rld, file = output_rld)

# === RESULTS === #
res <- results(dds, contrast = contrast, 
               independentFiltering = FALSE, cooksCutoff = Inf)

res <- lfcShrink(dds, contrast = contrast, res = res, type = "normal")
saveRDS(res, file = output_res)
