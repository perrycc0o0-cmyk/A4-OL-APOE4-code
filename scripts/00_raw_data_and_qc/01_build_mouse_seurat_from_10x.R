#!/usr/bin/env Rscript

# Build the mouse Seurat object from Cell Ranger filtered-feature matrices.
# Harmony is grouped by biological sample identity; APOE genotype is retained
# as metadata and is not supplied as a correction variable.

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(harmony)
})

source(file.path("R", "load_config.R"))

sample_key_file <- file.path(
  "metadata",
  "mouse_sample_key_and_doublet_qc.tsv"
)
sample_key <- read.delim(
  sample_key_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
sample_key <- sample_key[sample_key$internal_sample_id != "TOTAL", , drop = FALSE]

required_key_columns <- c("internal_sample_id", "public_sample_id", "apoe_group")
if (!all(required_key_columns %in% colnames(sample_key))) {
  stop(
    "Sample key is missing required columns: ",
    paste(setdiff(required_key_columns, colnames(sample_key)), collapse = ", ")
  )
}

qc_dir <- file.path(A4OL_MOUSE_QC_ROOT, "threshold_qc")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(A4OL_RESULTS_ROOT, recursive = TRUE, showWarnings = FALSE)

thresholds <- list(
  nFeature_min = 700,
  nFeature_max = 5000,
  nCount_min = 500,
  nCount_max = 30000,
  percent_mt_max = 5,
  log10GenesPerUMI_min = 0.75
)

read_and_filter_sample <- function(i) {
  internal_id <- sample_key$internal_sample_id[i]
  public_id <- sample_key$public_sample_id[i]
  genotype <- sample_key$apoe_group[i]
  input_h5 <- file.path(
    A4OL_MOUSE_RAW_ROOT,
    paste0(internal_id, "_outs"),
    "filtered_feature_bc_matrix.h5"
  )

  if (!file.exists(input_h5)) stop("Input matrix not found: ", input_h5)
  counts <- Seurat::Read10X_h5(input_h5)
  if (is.list(counts)) {
    if (!"Gene Expression" %in% names(counts)) {
      stop("Cannot identify the Gene Expression matrix in: ", input_h5)
    }
    counts <- counts[["Gene Expression"]]
  }

  obj <- Seurat::CreateSeuratObject(
    counts = counts,
    project = internal_id,
    min.cells = 3,
    min.features = 300
  )
  obj$sample <- internal_id
  obj$mouse_id_public <- public_id
  obj$genotype <- genotype
  obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^mt-")
  obj[["percent.ribo"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^Rp[sl]")
  obj$log10GenesPerUMI <- log10(obj$nFeature_RNA) / log10(obj$nCount_RNA)

  keep <-
    obj$nFeature_RNA > thresholds$nFeature_min &
    obj$nFeature_RNA < thresholds$nFeature_max &
    obj$nCount_RNA > thresholds$nCount_min &
    obj$nCount_RNA < thresholds$nCount_max &
    obj$percent.mt < thresholds$percent_mt_max &
    obj$log10GenesPerUMI > thresholds$log10GenesPerUMI_min

  stats <- data.frame(
    internal_sample_id = internal_id,
    public_sample_id = public_id,
    apoe_group = genotype,
    cells_after_CreateSeuratObject = ncol(obj),
    cells_after_threshold_qc = sum(keep),
    cells_removed_by_threshold_qc = sum(!keep),
    stringsAsFactors = FALSE
  )

  list(object = subset(obj, cells = colnames(obj)[keep]), stats = stats)
}

filtered <- lapply(seq_len(nrow(sample_key)), read_and_filter_sample)
names(filtered) <- sample_key$internal_sample_id

qc_summary <- do.call(rbind, lapply(filtered, `[[`, "stats"))
write.table(
  qc_summary,
  file.path(qc_dir, "threshold_qc_cell_counts_by_sample.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  data.frame(
    threshold = names(thresholds),
    value = unlist(thresholds, use.names = FALSE)
  ),
  file.path(qc_dir, "thresholds.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

objects <- lapply(filtered, `[[`, "object")
obj <- merge(
  x = objects[[1]],
  y = objects[-1],
  add.cell.ids = names(objects),
  merge.data = FALSE
)

obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(
  obj,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)
obj <- ScaleData(
  obj,
  features = VariableFeatures(obj),
  vars.to.regress = c("nCount_RNA", "percent.mt"),
  verbose = FALSE
)
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 50, verbose = FALSE)

n_pcs <- min(30L, ncol(Embeddings(obj, "pca")))
obj <- harmony::RunHarmony(
  object = obj,
  group.by.vars = "sample",
  reduction.use = "pca",
  dims.use = seq_len(n_pcs),
  project.dim = FALSE,
  verbose = TRUE
)

obj <- FindNeighbors(obj, reduction = "harmony", dims = seq_len(n_pcs), verbose = FALSE)
for (resolution in c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0)) {
  obj <- FindClusters(
    obj,
    resolution = resolution,
    random.seed = 42,
    verbose = FALSE
  )
}
obj$seurat_clusters <- obj$RNA_snn_res.0.4
obj <- RunUMAP(
  obj,
  reduction = "harmony",
  dims = seq_len(n_pcs),
  reduction.name = "umap",
  reduction.key = "UMAP_",
  n.neighbors = 30L,
  min.dist = 0.3,
  metric = "cosine",
  seed.use = 42,
  verbose = FALSE
)

expected_post_threshold_cells <- 120057L
if (ncol(obj) != expected_post_threshold_cells) {
  warning(
    "Post-threshold cell count is ", ncol(obj),
    "; the archived analysis contained ", expected_post_threshold_cells,
    ". Check raw-data release and package versions."
  )
}

saveRDS(
  obj,
  file.path(A4OL_RESULTS_ROOT, "seurat_final.rds"),
  compress = TRUE
)
write.csv(
  obj@meta.data,
  file.path(qc_dir, "post_threshold_cell_metadata.csv"),
  row.names = TRUE
)
capture.output(sessionInfo(), file = file.path(qc_dir, "sessionInfo.txt"))

message("Saved post-threshold object: ", file.path(A4OL_RESULTS_ROOT, "seurat_final.rds"))
message("Cells after threshold QC: ", ncol(obj))
