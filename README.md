# aml-stress-response-transcriptomics
[![DOI](https://zenodo.org/badge/1102508711.svg)](https://doi.org/10.5281/zenodo.17690695)
# run_all.R
# Master pipeline runner for AML chromatin regulator analysis

message("=== AML chromatin regulator workflow: start ===")

# Optionally set project root
# setwd("~/path/to/aml-chromatin-regulators")

source("scripts/01_download_and_import_data.R")
source("scripts/02_preprocessing_and_normalization.R")
source("scripts/03_deseq2.R")
source("scripts/04_enrichment.R")
source("scripts/05_network_analysis.R")

message("=== AML chromatin regulator workflow: complete ===")
# scripts/01_download_and_import_data.R
# Download and import GEO / TCGA data

suppressPackageStartupMessages({
  library(GEOquery)
  library(SummarizedExperiment)
  library(dplyr)
  library(readr)
})

dir.create("data", showWarnings = FALSE)
dir.create("data/raw", showWarnings = FALSE)
dir.create("metadata", showWarnings = FALSE)

message("Downloading GEO datasets (if not already present)...")

get_geo_counts <- function(gse_id) {
  message("Fetching ", gse_id, " ...")
  gse <- getGEO(gse_id, GSEMatrix = TRUE)
  gse <- gse[[1]]

  expr_mat <- exprs(gse)
  meta     <- pData(gse)

  list(expr = expr_mat, meta = meta)
}

# Discovery datasets (raw-like matrices from GEO)
gse247301 <- get_geo_counts("GSE247301")
gse251728 <- get_geo_counts("GSE251728")
gse282105 <- get_geo_counts("GSE282105")

counts_GSE247301 <- as.matrix(gse247301$expr)
meta_GSE247301   <- gse247301$meta

counts_GSE251728 <- as.matrix(gse251728$expr)
meta_GSE251728   <- gse251728$meta

counts_GSE282105 <- as.matrix(gse282105$expr)
meta_GSE282105   <- gse282105$meta

# Save raw matrices and metadata
write_csv(
  as.data.frame(counts_GSE247301) |>
    tibble::rownames_to_column("gene_id"),
  "data/raw/GSE247301_counts.csv"
)
write_csv(
  as.data.frame(counts_GSE251728) |>
    tibble::rownames_to_column("gene_id"),
  "data/raw/GSE251728_counts.csv"
)
write_csv(
  as.data.frame(counts_GSE282105) |>
    tibble::rownames_to_column("gene_id"),
  "data/raw/GSE282105_counts.csv"
)

write_csv(meta_GSE247301, "metadata/meta_GSE247301.csv")
write_csv(meta_GSE251728, "metadata/meta_GSE251728.csv")
write_csv(meta_GSE282105, "metadata/meta_GSE282105.csv")

message("GEO raw matrices and metadata saved in data/raw and metadata/.")

# Placeholder: TCGA-LAML and GSE229032
# In practice you will implement your own download or manual placement,
# then read them in here as matrices / data.frames.

message("Expecting TCGA-LAML and GSE229032 expression files to be present.")
# scripts/02_preprocessing_and_normalization.R
# Filter + normalize discovery datasets from raw counts

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(readr)
})

dir.create("data/processed", showWarnings = FALSE)

message("Preprocessing and normalization of GSE247301, GSE251728, GSE282105...")

normalize_dataset <- function(counts_mat, meta_df, dataset_id) {
  message("Normalizing ", dataset_id, " ...")

  # Ensure sample order matches
  meta_df <- meta_df[colnames(counts_mat), , drop = FALSE]

  # This assumes meta has a column 'condition'; adjust if needed
  if (!"condition" %in% colnames(meta_df)) {
    stop("Metadata for ", dataset_id, " must contain a 'condition' column.")
  }

  meta_df$condition <- factor(meta_df$condition)

  dds <- DESeqDataSetFromMatrix(
    countData = counts_mat,
    colData   = meta_df,
    design    = ~ condition
  )

  # Filter low counts
  keep <- rowSums(counts(dds) >= 10) >= 2
  dds  <- dds[keep, ]

  # Estimate size factors and dispersions
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)

  # Variance stabilizing transformation
  vsd <- vst(dds, blind = TRUE)

  norm_counts <- assay(vsd)

  out_file <- file.path("data/processed", paste0(dataset_id, "_norm_vst.csv"))
  write_csv(
    as.data.frame(norm_counts) |>
      tibble::rownames_to_column("gene_id"),
    out_file
  )

  message("Saved normalized VST counts for ", dataset_id, " -> ", out_file)

  invisible(list(dds = dds, vsd = vsd, norm_counts = norm_counts))
}

# Use objects created in 01_download_and_import_data.R
norm_GSE247301 <- normalize_dataset(counts_GSE247301, meta_GSE247301, "GSE247301")
norm_GSE251728 <- normalize_dataset(counts_GSE251728, meta_GSE251728, "GSE251728")
norm_GSE282105 <- normalize_dataset(counts_GSE282105, meta_GSE282105, "GSE282105")

saveRDS(norm_GSE247301$dds, "data/processed/dds_GSE247301.rds")
saveRDS(norm_GSE251728$dds, "data/processed/dds_GSE251728.rds")
saveRDS(norm_GSE282105$dds, "data/processed/dds_GSE282105.rds")

message("Preprocessing and normalization complete.")
# scripts/03_deseq2.R
# Differential expression analysis using DESeq2

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(readr)
})

dir.create("results", showWarnings = FALSE)
dir.create("results/DEG_tables", showWarnings = FALSE)

message("Running DESeq2 for differential expression...")

run_deseq <- function(dds_file, dataset_id) {
  message("DESeq2 for ", dataset_id, " ...")

  dds <- readRDS(dds_file)
  dds <- DESeq(dds)

  res <- results(dds)
  # Assumes coef 2 is the condition effect; adjust if more complex design
  res <- lfcShrink(dds, coef = 2, type = "apeglm", res = res)

  res_df <- as.data.frame(res) |>
    tibble::rownames_to_column("gene_id") |>
    arrange(padj)

  out_file <- file.path("results/DEG_tables", paste0("DEG_", dataset_id, ".csv"))
  write_csv(res_df, out_file)

  message("DE results saved -> ", out_file)

  invisible(res_df)
}

res_GSE247301 <- run_deseq("data/processed/dds_GSE247301.rds", "GSE247301")
res_GSE251728 <- run_deseq("data/processed/dds_GSE251728.rds", "GSE251728")
res_GSE282105 <- run_deseq("data/processed/dds_GSE282105.rds", "GSE282105")

message("DESeq2 differential expression completed for all discovery datasets.")
# scripts/04_enrichment.R
# Functional enrichment for DEGs and 73-gene module

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
  library(readr)
})

dir.create("results/enrichment", showWarnings = FALSE)

message("Running enrichment analysis...")

load_deg <- function(dataset_id, padj_cutoff = 0.05, lfc_cutoff = 1) {
  deg_file <- file.path("results/DEG_tables", paste0("DEG_", dataset_id, ".csv"))
  deg <- read_csv(deg_file, show_col_types = FALSE)

  sig <- deg |>
    filter(!is.na(padj), padj < padj_cutoff, abs(log2FoldChange) >= lfc_cutoff)

  sig$gene_id
}

deg_247301_genes <- load_deg("GSE247301")
deg_251728_genes <- load_deg("GSE251728")
deg_282105_genes <- load_deg("GSE282105")

# Example: GO BP enrichment for one dataset (you can loop or combine)
ego_247301 <- enrichGO(
  gene          = deg_247301_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

ego_file <- "results/enrichment/GO_BP_GSE247301.csv"
write_csv(as.data.frame(ego_247301), ego_file)
message("Saved GO BP enrichment for GSE247301 -> ", ego_file)

# Repeat similarly for KEGG, 73-gene module, etc. as needed.
message("Enrichment analysis script finished (extend for more datasets/modules).")
# scripts/05_network_analysis.R
# PPI network analysis and hub extraction

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

dir.create("results/ppi", showWarnings = FALSE)

message("Running PPI network analysis...")

ppi_file <- "73 ppi files .csv"  # or move into data/processed
ppi <- read_csv(ppi_file, show_col_types = FALSE)

# Basic degree-based hub list (EdgeCount)
degree_list <- ppi |>
  select(name, BetweennessCentrality, ClosenessCentrality, EdgeCount) |>
  arrange(desc(EdgeCount)) |>
  rename(
    degree = EdgeCount
  ) |>
  mutate(NumberofDirectedges = degree)

out_degree <- "results/ppi/73PPI_degree_ranked.csv"
write_csv(degree_list, out_degree)
message("Saved degree-ranked PPI hubs -> ", out_degree)

# Force-include CDKN1A, PHGDH, ALDH1L2 in top hubs
extra_hubs <- c("CDKN1A", "PHGDH", "ALDH1L2")

hub_rows <- degree_list |>
  filter(name %in% extra_hubs)

top10 <- degree_list |>
  slice(1:10)

combined <- bind_rows(top10, hub_rows) |>
  distinct(name, .keep_all = TRUE) |>
  arrange(desc(degree))

out_combined <- "results/ppi/73PPI_top10_plus_CDKN1A_PHGDH_ALDH1L2.csv"
write_csv(combined, out_combined)
message("Saved combined top10 + key hubs -> ", out_combined)

message("Network analysis script complete (STRING/Cytoscape steps documented in README).")

message("Expecting TCGA-LAML and GSE229032 expression files to be present in data/raw/")

# scripts/07_validation_GSE229032.R
suppressPackageStartupMessages({
  library(dplyr)
  library(pheatmap)
  library(ggplot2)
  library(ggrepel)
})

# assume: gse229032_counts, gse229032_meta, genes_73 are already in memory

aml_samples <- gse229032_meta |>
  filter(cell_type == "AML") |>      # example column
  pull(sample_id)

expr_aml <- gse229032_counts[, aml_samples]

expr_73 <- expr_aml[rownames(expr_aml) %in% genes_73, ]

expr_log2 <- log2(expr_73 + 1)

pheatmap(expr_log2, scale = "row", show_rownames = FALSE)

# PCA
pca <- prcomp(t(expr_log2), scale. = TRUE)
pca_df <- as.data.frame(pca$x) |>
  mutate(sample = colnames(expr_log2))

ggplot(pca_df, aes(PC1, PC2, label = sample)) +
  geom_point() +
  geom_text_repel()
  # scripts/06_integrative_signature_analysis.R
suppressPackageStartupMessages({
  library(GSVA)
  library(dplyr)
})

# assume: tcga_counts (genes x samples), tcga_meta, genes_73

gs <- gsva(
  as.matrix(tcga_counts),
  list(sig73 = genes_73),
  method = "gsva"
)

gsva_scores <- gs["sig73", ]
tcga_meta$GSVA_73 <- gsva_scores[tcga_meta$sample_id]

tcga_meta <- tcga_meta |>
  mutate(
    GSVA_73_group = ifelse(GSVA_73 >= median(GSVA_73, na.rm = TRUE), "High", "Low")
  )

saveRDS(tcga_meta, "data/processed/tcga_meta_with_GSVA.rds")
suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(dplyr)
  library(readr)
})

tcga_meta <- readRDS("data/processed/tcga_meta_with_GSVA.rds")

# K–M for GSVA
fit <- survfit(Surv(OS_time, OS_event) ~ GSVA_73_group, data = tcga_meta)

ggsurvplot(
  fit,
  data = tcga_meta,
  pval = TRUE,
  risk.table = TRUE
)

# Cox model
cox <- coxph(Surv(OS_time, OS_event) ~ GSVA_73_group, data = tcga_meta)
summary(cox)
TCGA-LAML and GSE229032 import is documented in:
`scripts/01_download_and_import_data.R`
(see the placeholder section for expected file formats and paths).

External validation of the 73-gene signature in GSE229032 (Methods 2.7) is implemented in:
`scripts/07_validation_GSE229032.R`.

GSVA scoring and survival analysis in TCGA-LAML (Methods 2.8) are implemented in:
`scripts/06_integrative_signature_analysis.R` and `scripts/08_figures.R`.
source("scripts/06_integrative_signature_analysis.R")
source("scripts/07_validation_GSE229032.R")
source("scripts/08_figures.R")
message("=== AML chromatin regulator workflow: complete ===")
