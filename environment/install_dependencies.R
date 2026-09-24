#!/usr/bin/env Rscript

snapshot_date <- "2026-02-28"
options(
  repos = c(
    CRAN = paste0("https://packagemanager.posit.co/cran/", snapshot_date)
  )
)

cran_packages <- c(
  "babelgene", "BiocManager", "cowplot", "data.table", "doParallel",
  "dplyr", "ggplot2",
  "ggrepel", "harmony", "hdf5r", "lmerTest", "patchwork", "pheatmap",
  "purrr", "readr", "remotes", "RColorBrewer", "scales", "Seurat",
  "SeuratObject", "showtext", "stringr", "sysfonts", "tibble", "tidyr"
)

if (.Platform$OS.type != "windows") cran_packages <- c(cran_packages, "doMC")

missing_cran <- cran_packages[
  !vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_cran) > 0) install.packages(missing_cran)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(version = "3.22", ask = FALSE)

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
  remotes::install_github(
    "aertslab/SCENIC@7a74341745cecd3505310c6c5755cad456756cf9"
  )
}

requirements <- read.delim(
  file.path("environment", "package_requirements.tsv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
if (.Platform$OS.type == "windows") {
  requirements <- requirements[requirements$package != "doMC", , drop = FALSE]
}

installed_versions <- vapply(
  requirements$package,
  function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) return(NA_character_)
    as.character(utils::packageVersion(pkg))
  },
  character(1)
)
expected_versions <- requirements$tested_version
mismatch <- is.na(installed_versions) | installed_versions != expected_versions

if (any(mismatch)) {
  warning(
    "The following direct dependencies do not match the declared 2026-02-28 ",
    "tested versions:\n",
    paste(
      sprintf(
        "  %s: expected %s; installed %s",
        requirements$package[mismatch],
        expected_versions[mismatch],
        ifelse(is.na(installed_versions[mismatch]), "not installed", installed_versions[mismatch])
      ),
      collapse = "\n"
    ),
    call. = FALSE
  )
}

message(
  "Dependency installation complete for the 2026-02-28 CRAN snapshot and ",
  "Bioconductor 3.22. Run environment/capture_session_info.R next."
)
