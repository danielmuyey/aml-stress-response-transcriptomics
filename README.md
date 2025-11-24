# aml-stress-response-transcriptomics
[![DOI](https://zenodo.org/badge/1102508711.svg)](https://doi.org/10.5281/zenodo.17690695)
# run_all.R
# Master script to run the full AML stress-response transcriptomics workflow

message("=== AML stress-response transcriptomics workflow: start ===")

# Discovery, DE, enrichment, network
source("scripts/01_download_and_import_data.R")
source("scripts/02_preprocessing_and_normalization.R")
source("scripts/03_deseq2.R")
source("scripts/04_enrichment.R")
source("scripts/05_network_analysis.R")

# Integrative TCGA-LAML analysis, external validation, figures
source("scripts/06_integrative_signature_analysis.R")  # GSVA + survival setup (TCGA-LAML)
source("scripts/07_validation_GSE229032.R")           # GSE229032 AML cell line validation
source("scripts/08_figures.R")                        # Final figures (heatmaps, PCA, KM, etc.)

message("=== AML chromatin regulator workflow: complete ===")
source("run_all.R")
# scripts/06_integrative_signature_analysis.R
# GSVA of 73-gene signature and hub genes in TCGA-LAML + survival data preparation

suppressPackageStartupMessages({
  library(GSVA)
  library(dplyr)
  library(readr)
})

message(">>> [06] Integrative signature analysis: TCGA-LAML GSVA and survival setup")

dir.create("data/processed", showWarnings = FALSE)
dir.create("results/survival", showWarnings = FALSE)

# Assumes tcga_counts, tcga_meta, and genes_73 are already available
# (or load them here if you're not keeping them in memory)

# Example loading if you use CSVs:
# tcga_counts <- read_csv("data/raw/TCGA_LAML_counts.csv") |>
#   tibble::column_to_rownames("gene_id") |>
#   as.matrix()
# tcga_meta <- read_csv("data/raw/TCGA_LAML_meta.csv")
# genes_73 <- read_csv("data/processed/73_gene_signature.csv")$gene_symbol

common_genes <- intersect(rownames(tcga_counts), genes_73)
if (length(common_genes) < 10) {
  warning("Fewer than 10 genes from the 73-gene signature are present in TCGA-LAML; check gene IDs.")
}

tcga_counts_sub <- tcga_counts[common_genes, , drop = FALSE]

gsva_res <- gsva(
  as.matrix(tcga_counts_sub),
  list(sig73 = common_genes),
  method = "gsva"
)

gsva_scores <- as.numeric(gsva_res["sig73", ])
names(gsva_scores) <- colnames(tcga_counts_sub)

tcga_meta <- tcga_meta |>
  mutate(
    GSVA_73 = gsva_scores[match(sample_id, names(gsva_scores))],
    GSVA_73_group = ifelse(GSVA_73 >= median(GSVA_73, na.rm = TRUE), "High", "Low")
  )

saveRDS(gsva_res, "data/processed/tcga_gsva_73_signature.rds")
saveRDS(tcga_meta, "data/processed/tcga_meta_with_GSVA.rds")

write_csv(
  tcga_meta,
  "results/survival/tcga_meta_with_GSVA.csv"
)

message("<<< [06] Integrative signature analysis completed; outputs written to 'data/processed/' and 'results/survival/'.")
# scripts/07_validation_GSE229032.R
# External validation of the 73-gene signature in GSE229032 AML cell lines

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(pheatmap)
  library(ggplot2)
  library(ggrepel)
})

message(">>> [07] External validation using GSE229032 AML cell lines")

dir.create("results/validation", showWarnings = FALSE)
dir.create("figs", showWarnings = FALSE)

# Assumes gse229032_counts, gse229032_meta, and genes_73 exist
# or load from CSVs if not in memory:

# gse229032_counts <- read_csv("data/raw/GSE229032_counts.csv") |>
#   tibble::column_to_rownames("gene_id") |>
#   as.matrix()
# gse229032_meta <- read_csv("data/raw/GSE229032_meta.csv")

# Filter to AML cell lines (adjust column names/values to your metadata)
aml_samples <- gse229032_meta |>
  filter(disease == "AML") |>           # or your real column/label
  slice_head(n = 15) |>                 # 15 AML lines as per Methods 2.7
  pull(sample_id)

expr_aml <- gse229032_counts[, aml_samples, drop = FALSE]

common_genes <- intersect(rownames(expr_aml), genes_73)
expr_73 <- expr_aml[common_genes, , drop = FALSE]

expr_log2 <- log2(expr_73 + 1)

# Heatmap
pdf("figs/GSE229032_73gene_heatmap.pdf", width = 6, height = 8)
pheatmap(
  expr_log2,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  show_rownames = FALSE,
  main = "GSE229032: 73-gene signature across AML cell lines"
)
dev.off()

# PCA
pca <- prcomp(t(expr_log2), scale. = TRUE)
pca_df <- as.data.frame(pca$x) |>
  mutate(sample_id = colnames(expr_log2))

pdf("figs/GSE229032_73gene_PCA.pdf", width = 6, height = 5)
ggplot(pca_df, aes(PC1, PC2, label = sample_id)) +
  geom_point(size = 2) +
  geom_text_repel(size = 2.5) +
  theme_minimal() +
  labs(
    title = "PCA of 73-gene signature in GSE229032 AML cell lines",
    x = "PC1",
    y = "PC2"
  )
dev.off()

write_csv(
  as.data.frame(expr_log2),
  "results/validation/GSE229032_73gene_log2expr.csv"
)

message("<<< [07] GSE229032 validation completed; heatmap/PCA saved to 'figs/' and log2 expression to 'results/validation/'.")

# ===============================================
# Placeholder: TCGA-LAML and GSE229032 datasets
# ===============================================

# In practice, you can either:
# 1) Programmatically download TCGA-LAML expression data (e.g., using TCGAbiolinks)
# 2) Place pre-processed CSV files manually in the data/raw/ directory

# Expected file locations:
# data/raw/TCGA_LAML_counts.csv      -> gene x sample count matrix
# data/raw/TCGA_LAML_meta.csv        -> metadata (sample info, condition, etc.)
# data/raw/GSE229032_counts.csv      -> gene x sample count matrix
# data/raw/GSE229032_meta.csv        -> metadata

# Example reading (once files are placed/downloaded):
# tcga_counts <- read_csv("data/raw/TCGA_LAML_counts.csv") |>
#   tibble::column_to_rownames("gene_id")
# tcga_meta   <- read_csv("data/raw/TCGA_LAML_meta.csv")
#
# gse229032_counts <- read_csv("data/raw/GSE229032_counts.csv") |>
#   tibble::column_to_rownames("gene_id")
# gse229032_meta   <- read_csv("data/raw/GSE229032_meta.csv")
## Reproducibility and reviewer notes

This repository is designed to make all analyses fully reproducible:

- **Scripts and pipeline:** All steps are implemented as modular R scripts in `scripts/` and orchestrated via `run_all.R`. Each script is documented in the README and in its header.
- **TCGA-LAML and GSE229032 placeholders:**
  - TCGA-LAML and GSE229032 data can be downloaded programmatically or supplied as CSVs in `data/raw/`.
  - Expected file names and example loading code are documented in `scripts/01_download_and_import_data.R`.
- **Pipeline diagram:** The end-to-end workflow from raw counts → normalization → DE → enrichment → network → validation → figures is summarized in the Mermaid diagram above.
- **Directory structure:** The `data/`, `metadata/`, `results/`, and `scripts/` directories are organized to match the analysis steps and figure/table outputs.
- **DOI and versioning:** This repository is archived on Zenodo and can be cited via the DOI badge at the top of this README, ensuring that the exact version used for the manuscript is preserved.

To reproduce the full analysis, from discovery datasets through external and translational validation, run:

```r
source("run_all.R")

