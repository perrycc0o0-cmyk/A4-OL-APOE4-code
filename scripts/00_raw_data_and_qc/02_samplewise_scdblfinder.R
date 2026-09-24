#!/usr/bin/env Rscript

# Detect doublets independently within each biological sample. The labelled
# object is retained, and only singlets are written to the downstream object.

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(scDblFinder)
  library(BiocParallel)
})

source(file.path("R", "load_config.R"))

input_rds <- file.path(A4OL_RESULTS_ROOT, "seurat_final.rds")
out_dir <- file.path(A4OL_MOUSE_QC_ROOT, "scDblFinder_by_sample")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_rds)) stop("Input object not found: ", input_rds)
obj <- readRDS(input_rds)
if (!"RNA" %in% names(obj@assays)) stop("RNA assay is absent from input object.")
if (!"sample" %in% colnames(obj@meta.data)) {
  stop("Metadata column 'sample' is required for sample-wise doublet detection.")
}
if (anyNA(obj$sample) || any(!nzchar(as.character(obj$sample)))) {
  stop("Metadata column 'sample' contains missing or empty values.")
}

DefaultAssay(obj) <- "RNA"
obj_work <- obj
if (exists(
  "JoinLayers",
  envir = asNamespace("SeuratObject"),
  mode = "function",
  inherits = FALSE
)) {
  obj_work <- tryCatch(
    SeuratObject::JoinLayers(obj_work, assay = "RNA"),
    error = function(e) obj_work
  )
}

sce <- Seurat::as.SingleCellExperiment(obj_work, assay = "RNA")
rm(obj_work)
SummarizedExperiment::colData(sce)$doublet_qc_sample <- as.character(obj$sample)

if (.Platform$OS.type == "unix" && A4OL_N_CORES > 1L) {
  bp <- BiocParallel::MulticoreParam(
    workers = A4OL_N_CORES,
    RNGseed = 20260820,
    progressbar = TRUE
  )
} else {
  # SerialParam is supported on Windows and avoids fork-specific behavior.
  bp <- BiocParallel::SerialParam(
    RNGseed = 20260820,
    progressbar = TRUE
  )
}

sce <- scDblFinder::scDblFinder(
  sce,
  samples = "doublet_qc_sample",
  clusters = TRUE,
  dbr = NULL,
  multiSampleMode = "split",
  BPPARAM = bp
)

sc_cols <- grep(
  "^scDblFinder\\.",
  colnames(SummarizedExperiment::colData(sce)),
  value = TRUE
)
if (!all(c("scDblFinder.score", "scDblFinder.class") %in% sc_cols)) {
  stop("scDblFinder did not return the required score and class columns.")
}
for (column in sc_cols) {
  value <- SummarizedExperiment::colData(sce)[[column]]
  names(value) <- colnames(sce)
  obj[[column]] <- value[colnames(obj)]
}
rm(sce)

classes <- as.character(obj$scDblFinder.class)
summary_by_sample <- do.call(
  rbind,
  lapply(split(seq_len(ncol(obj)), as.character(obj$sample)), function(index) {
    local_class <- classes[index]
    data.frame(
      sample = as.character(obj$sample[index[1]]),
      cells_before_doublet_qc = length(index),
      doublets_removed = sum(local_class == "doublet"),
      singlets_retained = sum(local_class == "singlet"),
      doublet_rate = mean(local_class == "doublet"),
      stringsAsFactors = FALSE
    )
  })
)

singlet_cells <- colnames(obj)[classes == "singlet"]
doublet_cells <- colnames(obj)[classes == "doublet"]
obj_singlet <- subset(obj, cells = singlet_cells)

write.table(
  summary_by_sample,
  file.path(out_dir, "doublet_summary_by_sample.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
writeLines(singlet_cells, file.path(out_dir, "singlet_cell_barcodes.txt"))
writeLines(doublet_cells, file.path(out_dir, "doublet_cell_barcodes_removed.txt"))
saveRDS(obj, file.path(out_dir, "seurat_final_with_scDblFinder_labels.rds"))
saveRDS(
  obj_singlet,
  file.path(out_dir, "seurat_final_singlet_only_scDblFinder_by_sample.rds")
)

# The annotation stage reads this singlet-only downstream object.
saveRDS(
  obj_singlet,
  file.path(A4OL_RESULTS_ROOT, "seurat_clustered.rds"),
  compress = TRUE
)
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))

expected <- c(total = 120057L, singlet = 108166L, doublet = 11891L)
observed <- c(
  total = ncol(obj),
  singlet = length(singlet_cells),
  doublet = length(doublet_cells)
)
if (!identical(unname(observed), unname(expected))) {
  warning(
    "Observed counts differ from the archived analysis: ",
    paste(names(observed), observed, sep = "=", collapse = ", ")
  )
}

message("Input cells: ", observed[["total"]])
message("Doublets removed: ", observed[["doublet"]])
message("Singlets retained for downstream analysis: ", observed[["singlet"]])
