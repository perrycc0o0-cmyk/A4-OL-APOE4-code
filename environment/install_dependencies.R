#!/usr/bin/env Rscript

options(repos = c(CRAN = "https://cloud.r-project.org"))

cran_packages <- c(
  "babelgene", "cowplot", "data.table", "doParallel", "dplyr", "ggplot2",
  "ggrepel", "harmony", "hdf5r", "lmerTest", "patchwork", "pheatmap",
  "purrr", "readr", "remotes", "RColorBrewer", "scales", "Seurat",
  "SeuratObject", "showtext", "stringr", "sysfonts", "tibble", "tidyr"
)

if (.Platform$OS.type != "windows") cran_packages <- c(cran_packages, "doMC")

missing_cran <- cran_packages[
  !vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_cran) > 0) install.packages(missing_cran)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c(
  "AnnotationDbi", "AUCell", "BiocParallel", "clusterProfiler", "enrichplot",
  "fgsea", "GENIE3", "org.Hs.eg.db", "org.Mm.eg.db", "RcisTarget",
  "scDblFinder", "SingleCellExperiment", "SummarizedExperiment", "UCell"
)
missing_bioc <- bioc_packages[
  !vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_bioc) > 0) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

if (!requireNamespace("SCENIC", quietly = TRUE)) {
  remotes::install_github("aertslab/SCENIC")
}

message("Dependency installation complete. Run environment/capture_session_info.R next.")
