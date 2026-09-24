#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# ============================================================
# STEP2_SCENIC_A4OL_FORMAL_HUMAN_MOUSE_V9.R
#
# Formal SCENIC run for A4-OLs:
#   1) Human A4-like OLs vs Other OLs
#   2) Mouse A4-OLs vs Other OLs
#
# This version incorporates fixes verified by the minimal test:
#   - v9 gene-based cisTarget database uses index column: motifs
#   - use motifAnnotations_hgnc_v9 / motifAnnotations_mgi_v9
#   - bypass SCENIC::geneFiltering compatibility issues
#   - save cellInfo as RDS for SCENIC heatmap step
#   - use doMC for AUCell / SCENIC multicore scoring
#   - keep EGR1/Egr1 in expression matrix when present
# ============================================================

options(stringsAsFactors = FALSE)
options(timeout = 300)
options(repos = c(CRAN = "https://cloud.r-project.org"))

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(data.table)
  library(SCENIC)
  library(RcisTarget)
  library(AUCell)
  library(GENIE3)
  library(doMC)
})

source(file.path("R", "load_config.R"))

HAS_GGPLOT2 <- requireNamespace("ggplot2", quietly = TRUE)

cat0 <- function(...) cat(..., "\n", sep = "")

# ============================================================
# 0. User parameters
# ============================================================
SEED <- 123
set.seed(SEED)

# Your machine has 32 cores. 24 is more stable than using all cores.
N_CORES <- 30

doMC::registerDoMC(cores = N_CORES)

# Formal but still stable settings.
# If you want a faster pilot, use 1000 and 8000.
# If you want a larger final run, use 3000 and 12000/Inf.
MAX_CELLS_PER_GROUP <- 1500
MAX_GENES <- 10000

# Sparse gene filter. EGR1/Egr1 will be force-kept if present and detected.
MIN_CELLS_FACTOR <- 0.01
MIN_COUNTS_FACTOR <- 0.03
MIN_CELLS_ABSOLUTE <- 3
MIN_COUNTS_ABSOLUTE <- 3

# Set SCENIC_CLEAN_OUTPUT=true only when an existing species output directory
# should be removed before a fresh run. The default preserves existing results.
CLEAN_SPECIES_OUTPUT_DIR <- tolower(
  Sys.getenv("SCENIC_CLEAN_OUTPUT", unset = "false")
) %in% c("true", "1", "yes")

BASE_DIR <- A4OL_SERVER_ROOT

OUT_ROOT <- file.path(A4OL_FIGURE_ROOT, "SCENIC_A4OL_FORMAL_V9")

# Input RDS files
HUMAN_RDS <- file.path(
  A4OL_FIGURE_ROOT,
  "Mapping",
  "res",
  "res_0.4",
  "oligodendrocytes_res0.4.rds"
)

MOUSE_RDS <- file.path(
  A4OL_FIGURE_ROOT,
  "olig-c",
  "seurat_oligodendrocytes_subset_annotated_A4_OLs.rds"
)

# cisTarget databases
DB_DIR_HUMAN <- file.path(A4OL_CISTARGET_ROOT, "hg38")
DB_DIR_MOUSE <- file.path(A4OL_CISTARGET_ROOT, "mm10")

DBS_HUMAN <- c(
  "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather",
  "hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
)

DBS_MOUSE <- c(
  "mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather",
  "mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
)

# Force keep these TF genes if present in the expression object.
FORCE_KEEP_HUMAN <- c("EGR1")
FORCE_KEEP_MOUSE <- c("Egr1")

# ============================================================
# 1. Compatibility patch for v9 gene-based cisTarget databases
# ============================================================
patch_scenic_for_v9_gene_based_db <- function(db_index_col = "motifs") {
  ns_scenic <- asNamespace("SCENIC")
  ns_rct <- asNamespace("RcisTarget")

  patched_getDbAnnotations <- function(scenicOptions) {
    dbAnnotFiles <- scenicOptions@settings$db_annotFiles

    if (!is.null(dbAnnotFiles)) {
      motifAnnotations <- NULL
      for (annotPath in dbAnnotFiles) {
        motifAnnot <- data.table::fread(annotPath)
        motifAnnot$annotationSource <- factor(motifAnnot$annotationSource)
        colnames(motifAnnot)[1] <- "motif"
        motifAnnotations <- rbind(motifAnnotations, motifAnnot)
      }
      return(motifAnnotations)
    }

    org <- SCENIC::getDatasetInfo(scenicOptions, "org")
    if (is.na(org) || !org %in% c("hgnc", "mgi", "dmel")) {
      stop("Organism not recognized: ", org)
    }

    motifAnnotCandidates <- switch(
      org,
      hgnc = c("motifAnnotations_hgnc_v9", "motifAnnotations_hgnc", "motifAnnotations"),
      mgi  = c("motifAnnotations_mgi_v9",  "motifAnnotations_mgi",  "motifAnnotations"),
      dmel = c("motifAnnotations_dmel_v9", "motifAnnotations_dmel", "motifAnnotations")
    )

    for (motifAnnotName in motifAnnotCandidates) {
      env <- new.env(parent = emptyenv())

      ok <- suppressWarnings(
        tryCatch({
          utils::data(list = motifAnnotName, package = "RcisTarget", envir = env)
          TRUE
        }, error = function(e) FALSE)
      )

      if (ok && exists(motifAnnotName, envir = env, inherits = FALSE)) {
        return(get(motifAnnotName, envir = env, inherits = FALSE))
      }

      if (ok && exists("motifAnnotations", envir = env, inherits = FALSE)) {
        return(get("motifAnnotations", envir = env, inherits = FALSE))
      }
    }

    stop(
      "Cannot load motif annotations from RcisTarget for org = ", org,
      ". Tried: ", paste(motifAnnotCandidates, collapse = ", ")
    )
  }

  patched_checkAnnots <- function(object, motifAnnot) {
    allFeaturesInAnnot <- unique(as.character(unlist(motifAnnot[, 1])))

    featuresWithAnnot <- lapply(SCENIC::getDatabases(object), function(dbFile) {
      rnks <- tryCatch({
        RcisTarget::getRowNames(dbFile, indexCol = db_index_col)
      }, error = function(e) {
        RcisTarget::getRowNames(dbFile, indexCol = NULL)
      })

      rnks <- unique(as.character(rnks))
      length(intersect(allFeaturesInAnnot, rnks)) /
        length(unique(c(allFeaturesInAnnot, rnks)))
    })

    return(featuresWithAnnot)
  }

  original_importRankings <- get("importRankings", envir = ns_rct)

  patched_importRankings <- function(dbFilePath,
                                     columns = NULL,
                                     indexCol = db_index_col,
                                     warnMissingColumns = TRUE,
                                     ...) {
    original_importRankings(
      dbFilePath,
      columns = columns,
      indexCol = db_index_col,
      warnMissingColumns = warnMissingColumns,
      ...
    )
  }

  for (fn in c("getDbAnnotations", "checkAnnots")) {
    if (exists(fn, envir = ns_scenic, inherits = FALSE)) {
      if (bindingIsLocked(fn, ns_scenic)) unlockBinding(fn, ns_scenic)
      assign(
        fn,
        if (fn == "getDbAnnotations") patched_getDbAnnotations else patched_checkAnnots,
        envir = ns_scenic
      )
      lockBinding(fn, ns_scenic)
    }
  }

  if (exists("importRankings", envir = ns_scenic, inherits = FALSE)) {
    if (bindingIsLocked("importRankings", ns_scenic)) unlockBinding("importRankings", ns_scenic)
    assign("importRankings", patched_importRankings, envir = ns_scenic)
    lockBinding("importRankings", ns_scenic)
  }

  if (bindingIsLocked("importRankings", ns_rct)) unlockBinding("importRankings", ns_rct)
  assign("importRankings", patched_importRankings, envir = ns_rct)
  lockBinding("importRankings", ns_rct)

  cat0("[OK] SCENIC v9 compatibility patch applied. dbIndexCol = ", db_index_col)
}

load_motif_annotation_aliases <- function() {
  cat0("\n========== Load motif annotations ==========")

  ok_hgnc_v9 <- tryCatch({
    data("motifAnnotations_hgnc_v9", package = "RcisTarget", envir = .GlobalEnv)
    TRUE
  }, error = function(e) FALSE)

  if (ok_hgnc_v9 && exists("motifAnnotations_hgnc_v9", envir = .GlobalEnv)) {
    assign("motifAnnotations_hgnc", get("motifAnnotations_hgnc_v9", envir = .GlobalEnv), envir = .GlobalEnv)
    cat0("[OK] motifAnnotations_hgnc_v9 loaded and aliased to motifAnnotations_hgnc.")
  } else {
    ok_hgnc_old <- tryCatch({
      data("motifAnnotations_hgnc", package = "RcisTarget", envir = .GlobalEnv)
      TRUE
    }, error = function(e) FALSE)
    if (!ok_hgnc_old) stop("Cannot load motifAnnotations_hgnc_v9 or motifAnnotations_hgnc.")
    cat0("[OK] motifAnnotations_hgnc loaded.")
  }

  ok_mgi_v9 <- tryCatch({
    data("motifAnnotations_mgi_v9", package = "RcisTarget", envir = .GlobalEnv)
    TRUE
  }, error = function(e) FALSE)

  if (ok_mgi_v9 && exists("motifAnnotations_mgi_v9", envir = .GlobalEnv)) {
    assign("motifAnnotations_mgi", get("motifAnnotations_mgi_v9", envir = .GlobalEnv), envir = .GlobalEnv)
    cat0("[OK] motifAnnotations_mgi_v9 loaded and aliased to motifAnnotations_mgi.")
  } else {
    ok_mgi_old <- tryCatch({
      data("motifAnnotations_mgi", package = "RcisTarget", envir = .GlobalEnv)
      TRUE
    }, error = function(e) FALSE)
    if (!ok_mgi_old) stop("Cannot load motifAnnotations_mgi_v9 or motifAnnotations_mgi.")
    cat0("[OK] motifAnnotations_mgi loaded.")
  }
}

cat0("\n========== Apply SCENIC v9 patch ==========")
patch_scenic_for_v9_gene_based_db("motifs")
load_motif_annotation_aliases()

# ============================================================
# 2. Helper functions
# ============================================================
check_db_files <- function(db_dir, dbs, label) {
  cat0("\n========== Check ", label, " cisTarget databases ==========")
  db_files <- file.path(db_dir, dbs)
  db_info <- data.frame(
    file = db_files,
    exists = file.exists(db_files),
    size_gb = round(file.size(db_files) / 1024^3, 3),
    stringsAsFactors = FALSE
  )
  print(db_info)

  if (any(!db_info$exists)) stop(label, ": some database files do not exist.")
  if (any(is.na(db_info$size_gb)) || any(db_info$size_gb <= 0)) {
    stop(label, ": some database files are empty or invalid.")
  }

  cat0("Test importRankings(indexCol = motifs) for ", label)
  test_db <- RcisTarget::importRankings(
    dbFilePath = db_files[1],
    columns = NULL,
    indexCol = "motifs"
  )
  print(class(test_db))
  rm(test_db)
  gc()
  invisible(db_info)
}

get_counts_matrix <- function(obj) {
  counts <- tryCatch({
    GetAssayData(obj, assay = DefaultAssay(obj), layer = "counts")
  }, error = function(e) {
    GetAssayData(obj, assay = DefaultAssay(obj), slot = "counts")
  })
  return(counts)
}

strip_regulon_tf <- function(x) {
  y <- gsub("_extended", "", x)
  y <- gsub("\\s*\\(.*\\)$", "", y)
  y <- trimws(y)
  y
}

safe_wilcox <- function(x, group) {
  ok <- is.finite(x) & !is.na(group)
  x <- x[ok]
  group <- group[ok]
  if (length(unique(group)) != 2) return(NA_real_)
  if (length(x[group == unique(group)[1]]) < 3 || length(x[group == unique(group)[2]]) < 3) return(NA_real_)
  p <- tryCatch({
    wilcox.test(x ~ group)$p.value
  }, error = function(e) NA_real_)
  p
}

make_human_a4_status <- function(obj) {
  meta <- obj@meta.data
  if ("olig_clusters" %in% colnames(meta)) {
    obj$olig_clusters <- as.character(meta$olig_clusters)
  } else if ("RNA_snn_res.0.4" %in% colnames(meta)) {
    obj$olig_clusters <- as.character(meta$RNA_snn_res.0.4)
  } else if ("seurat_clusters" %in% colnames(meta)) {
    obj$olig_clusters <- as.character(meta$seurat_clusters)
  } else {
    stop("Human object: cannot find olig_clusters, RNA_snn_res.0.4, or seurat_clusters.")
  }

  # marker-0.4.R: clusters 2 + 12 are human A4-like OLs
  obj$A4_status <- ifelse(obj$olig_clusters %in% c("2", "12"), "A4-OLs", "Other OLs")
  obj$A4_status <- factor(obj$A4_status, levels = c("Other OLs", "A4-OLs"))
  obj
}

make_mouse_a4_status <- function(obj) {
  meta <- obj@meta.data

  candidate_cols <- c(
    "A4_status",
    "oligo_annotated_cluster",
    "olig_annotated_cluster",
    "annotated_cluster",
    "cluster_annotation",
    "manual_cluster",
    "cell_type_manual",
    "cell_type",
    "olig_clusters",
    "seurat_clusters"
  )

  found <- NULL
  for (cc in candidate_cols) {
    if (cc %in% colnames(meta)) {
      vals <- as.character(meta[[cc]])
      if (any(vals == "A4-OLs", na.rm = TRUE)) {
        found <- cc
        obj$A4_status <- ifelse(vals == "A4-OLs", "A4-OLs", "Other OLs")
        break
      }
    }
  }

  if (is.null(found)) {
    # Fallback based on previous olig_annotation-c.R definition: mouse olig clusters 1 + 4 = A4-OLs
    cluster_col <- NULL
    for (cc in c("olig_clusters", "seurat_clusters", "RNA_snn_res.0.4", "RNA_snn_res.0.5")) {
      if (cc %in% colnames(meta)) {
        cluster_col <- cc
        break
      }
    }
    if (is.null(cluster_col)) {
      stop("Mouse object: cannot detect A4 annotation or cluster column for fallback.")
    }
    vals <- as.character(meta[[cluster_col]])
    obj$A4_status <- ifelse(vals %in% c("1", "4"), "A4-OLs", "Other OLs")
    found <- paste0(cluster_col, " fallback: 1+4")
  }

  obj$A4_status <- factor(obj$A4_status, levels = c("Other OLs", "A4-OLs"))
  attr(obj$A4_status, "source_column") <- found
  obj
}

prepare_expr_for_scenic <- function(obj,
                                    label,
                                    out_dirs,
                                    force_keep_genes,
                                    max_cells_per_group,
                                    max_genes) {
  cat0("\n========== ", label, ": prepare cells and expression ==========")

  if ("RNA" %in% Assays(obj)) DefaultAssay(obj) <- "RNA"

  meta <- obj@meta.data
  if (!"A4_status" %in% colnames(meta)) stop(label, ": A4_status was not created.")

  cat0(label, " A4 status table before sampling:")
  print(table(meta$A4_status, useNA = "ifany"))

  cells_other <- rownames(meta)[meta$A4_status == "Other OLs"]
  cells_a4 <- rownames(meta)[meta$A4_status == "A4-OLs"]

  if (length(cells_other) < max_cells_per_group) {
    stop(label, ": Other OLs fewer than MAX_CELLS_PER_GROUP. n = ", length(cells_other))
  }
  if (length(cells_a4) < max_cells_per_group) {
    stop(label, ": A4-OLs fewer than MAX_CELLS_PER_GROUP. n = ", length(cells_a4))
  }

  cells_use <- c(
    sample(cells_other, max_cells_per_group),
    sample(cells_a4, max_cells_per_group)
  )

  obj_use <- subset(obj, cells = cells_use)
  meta_use <- obj_use@meta.data

  cat0(label, " cells used:")
  print(table(meta_use$A4_status))

  counts <- get_counts_matrix(obj_use)
  cat0(label, " raw counts dim:")
  print(dim(counts))

  min_cells <- max(MIN_CELLS_ABSOLUTE, ceiling(ncol(counts) * MIN_CELLS_FACTOR))
  min_counts <- max(MIN_COUNTS_ABSOLUTE, ceiling(ncol(counts) * MIN_COUNTS_FACTOR))

  gene_ncells <- Matrix::rowSums(counts > 0)
  gene_counts <- Matrix::rowSums(counts)

  sparse_keep <- gene_ncells >= min_cells & gene_counts >= min_counts
  sparse_genes <- rownames(counts)[sparse_keep]

  force_present <- intersect(force_keep_genes, rownames(counts))
  force_detected <- force_present[gene_counts[force_present] > 0]

  genes_after_sparse <- union(sparse_genes, force_detected)

  counts_f <- counts[genes_after_sparse, , drop = FALSE]

  cat0(label, " genes after sparse filter plus force-kept detected genes: ", nrow(counts_f))

  if (is.finite(max_genes) && nrow(counts_f) > max_genes) {
    gene_means <- Matrix::rowMeans(counts_f)
    gene_sq_means <- Matrix::rowMeans(counts_f ^ 2)
    gene_vars <- gene_sq_means - gene_means ^ 2

    top_genes <- names(sort(gene_vars, decreasing = TRUE))[seq_len(max_genes)]
    final_genes <- union(top_genes, force_detected)
    counts_f <- counts_f[final_genes, , drop = FALSE]
  }

  cat0(label, " genes after MAX_GENES limit plus force-kept detected genes: ", nrow(counts_f))

  exprMat <- as.matrix(log2(counts_f + 1))
  storage.mode(exprMat) <- "numeric"

  cat0(label, " exprMat dim:")
  print(dim(exprMat))

  force_report <- data.frame(
    species = label,
    gene = force_keep_genes,
    present_in_object = force_keep_genes %in% rownames(counts),
    detected_in_sampled_cells = force_keep_genes %in% force_detected,
    kept_in_exprMat = force_keep_genes %in% rownames(exprMat),
    stringsAsFactors = FALSE
  )

  write.csv(
    force_report,
    file = file.path(out_dirs$tables, paste0(label, "_force_keep_gene_report.csv")),
    row.names = FALSE
  )

  cellInfo <- data.frame(
    cell = rownames(meta_use),
    A4_status = as.character(meta_use$A4_status),
    stringsAsFactors = FALSE,
    row.names = rownames(meta_use)
  )

  # Add useful metadata if present.
  extra_cols <- intersect(
    c("sample", "orig.ident", "APOE", "apoe_group", "apoe_genotype", "diagnosis", "type", "dataset", "olig_clusters", "seurat_clusters"),
    colnames(meta_use)
  )
  if (length(extra_cols) > 0) {
    cellInfo <- cbind(cellInfo, meta_use[rownames(cellInfo), extra_cols, drop = FALSE])
  }

  write.csv(
    cellInfo,
    file = file.path(out_dirs$tables, paste0(label, "_cellInfo.csv")),
    row.names = FALSE,
    quote = FALSE
  )

  saveRDS(
    exprMat,
    file = file.path(out_dirs$rds, paste0(label, "_exprMat_log2_counts.rds")),
    compress = TRUE
  )

  saveRDS(
    cellInfo,
    file = file.path(out_dirs$rds, paste0(label, "_cellInfo.rds")),
    compress = TRUE
  )

  gene_filter_summary <- data.frame(
    species = label,
    n_cells_total = ncol(obj),
    n_cells_used = ncol(exprMat),
    n_other_used = sum(cellInfo$A4_status == "Other OLs"),
    n_a4_used = sum(cellInfo$A4_status == "A4-OLs"),
    n_genes_raw = nrow(counts),
    n_genes_sparse_plus_force = length(genes_after_sparse),
    n_genes_final = nrow(exprMat),
    min_cells = min_cells,
    min_counts = min_counts,
    max_genes = max_genes,
    stringsAsFactors = FALSE
  )

  write.csv(
    gene_filter_summary,
    file = file.path(out_dirs$tables, paste0(label, "_SCENIC_input_summary.csv")),
    row.names = FALSE
  )

  rm(counts, counts_f, obj_use)
  gc()

  list(exprMat = exprMat, cellInfo = cellInfo, force_report = force_report)
}

calc_regulon_stats <- function(auc_mat, cellInfo, label, out_dirs) {
  common_cells <- intersect(colnames(auc_mat), rownames(cellInfo))
  if (length(common_cells) == 0) stop(label, ": no common cells between AUC matrix and cellInfo.")

  auc_mat <- auc_mat[, common_cells, drop = FALSE]
  cellInfo <- cellInfo[common_cells, , drop = FALSE]

  group <- factor(cellInfo$A4_status, levels = c("Other OLs", "A4-OLs"))
  regs <- rownames(auc_mat)

  stats_list <- lapply(regs, function(reg) {
    x <- as.numeric(auc_mat[reg, ])
    other <- x[group == "Other OLs"]
    a4 <- x[group == "A4-OLs"]
    p <- safe_wilcox(x, group)
    data.frame(
      species = label,
      regulon = reg,
      TF = strip_regulon_tf(reg),
      TF_upper = toupper(strip_regulon_tf(reg)),
      n_other = length(other),
      n_a4 = length(a4),
      mean_other = mean(other, na.rm = TRUE),
      mean_a4 = mean(a4, na.rm = TRUE),
      median_other = median(other, na.rm = TRUE),
      median_a4 = median(a4, na.rm = TRUE),
      delta_mean_A4_minus_other = mean(a4, na.rm = TRUE) - mean(other, na.rm = TRUE),
      wilcox_p = p,
      stringsAsFactors = FALSE
    )
  })

  stats <- bind_rows(stats_list)
  stats$wilcox_fdr <- p.adjust(stats$wilcox_p, method = "BH")
  stats <- stats %>% arrange(desc(delta_mean_A4_minus_other), wilcox_fdr)

  write.csv(
    stats,
    file = file.path(out_dirs$tables, paste0(label, "_all_regulon_AUC_A4_vs_Other_stats.csv")),
    row.names = FALSE
  )

  egr1_stats <- stats %>% filter(TF_upper == "EGR1")
  write.csv(
    egr1_stats,
    file = file.path(out_dirs$tables, paste0(label, "_EGR1_regulon_AUC_stats.csv")),
    row.names = FALSE
  )

  if (nrow(egr1_stats) > 0) {
    egr1_reg <- egr1_stats$regulon[1]
    egr1_df <- data.frame(
      cell = common_cells,
      A4_status = as.character(group),
      EGR1_regulon_AUC = as.numeric(auc_mat[egr1_reg, common_cells]),
      regulon_name = egr1_reg,
      stringsAsFactors = FALSE
    )

    write.csv(
      egr1_df,
      file = file.path(out_dirs$tables, paste0(label, "_EGR1_regulon_AUC_per_cell.csv")),
      row.names = FALSE
    )

    if (HAS_GGPLOT2) {
      p <- ggplot2::ggplot(egr1_df, ggplot2::aes(x = A4_status, y = EGR1_regulon_AUC, fill = A4_status)) +
        ggplot2::geom_violin(trim = FALSE, alpha = 0.25, linewidth = 0.2) +
        ggplot2::geom_boxplot(width = 0.22, outlier.shape = NA, alpha = 0.75, linewidth = 0.35) +
        ggplot2::geom_jitter(width = 0.15, size = 0.35, alpha = 0.35) +
        ggplot2::theme_classic(base_size = 12) +
        ggplot2::labs(
          title = paste0(label, " ", egr1_reg, " activity"),
          x = NULL,
          y = "Regulon AUC"
        ) +
        ggplot2::theme(legend.position = "none")

      ggplot2::ggsave(
        filename = file.path(out_dirs$plots, paste0(label, "_EGR1_regulon_AUC_A4_vs_Other.png")),
        plot = p,
        width = 4.2,
        height = 4.2,
        dpi = 300
      )
      ggplot2::ggsave(
        filename = file.path(out_dirs$plots, paste0(label, "_EGR1_regulon_AUC_A4_vs_Other.pdf")),
        plot = p,
        width = 4.2,
        height = 4.2
      )
    }
  } else {
    cat0(label, ": no EGR1/Egr1 regulon detected. This is a biological/threshold result, not a pipeline failure.")
  }

  if (HAS_GGPLOT2 && nrow(stats) > 0) {
    top_stats <- stats %>%
      filter(is.finite(delta_mean_A4_minus_other)) %>%
      arrange(desc(abs(delta_mean_A4_minus_other))) %>%
      head(25)

    top_stats$regulon <- factor(top_stats$regulon, levels = rev(top_stats$regulon))

    p2 <- ggplot2::ggplot(top_stats, ggplot2::aes(x = regulon, y = delta_mean_A4_minus_other)) +
      ggplot2::geom_col(width = 0.75) +
      ggplot2::coord_flip() +
      ggplot2::theme_classic(base_size = 10) +
      ggplot2::labs(
        title = paste0(label, " top regulon AUC difference"),
        x = NULL,
        y = "Mean AUC difference: A4-OLs - Other OLs"
      )

    ggplot2::ggsave(
      filename = file.path(out_dirs$plots, paste0(label, "_top25_regulon_delta_A4_minus_Other.png")),
      plot = p2,
      width = 7,
      height = 6,
      dpi = 300
    )
    ggplot2::ggsave(
      filename = file.path(out_dirs$plots, paste0(label, "_top25_regulon_delta_A4_minus_Other.pdf")),
      plot = p2,
      width = 7,
      height = 6
    )
  }

  stats
}

extract_auc_object <- function(scenicOptions, scenic_dir) {
  oldwd <- getwd()
  on.exit(setwd(oldwd), add = TRUE)
  setwd(scenic_dir)

  regulonAUC <- tryCatch({
    SCENIC::loadInt(scenicOptions, "aucell_regulonAUC")
  }, error = function(e) {
    auc_candidates <- c(
      file.path(scenic_dir, "int", "3.4_regulonAUC.Rds"),
      file.path(scenic_dir, "int", "3.4_regulonAUC.rds")
    )
    auc_file <- auc_candidates[file.exists(auc_candidates)][1]

    if (is.na(auc_file)) {
      all_rds <- list.files(file.path(scenic_dir, "int"), pattern = "\\.Rds$|\\.rds$", recursive = TRUE, full.names = TRUE)
      auc_file <- all_rds[grepl("regulonAUC|AUC", basename(all_rds), ignore.case = TRUE)][1]
    }

    if (is.na(auc_file) || !file.exists(auc_file)) {
      stop("Cannot find regulon AUC RDS under: ", file.path(scenic_dir, "int"))
    }
    readRDS(auc_file)
  })

  regulonAUC
}

run_scenic_one_species <- function(label,
                                   org,
                                   rds_file,
                                   db_dir,
                                   dbs,
                                   force_keep_genes,
                                   make_status_fun) {
  cat0("\n\n############################################################")
  cat0("Running formal SCENIC for: ", label)
  cat0("############################################################")

  species_root <- file.path(OUT_ROOT, label)
  if (CLEAN_SPECIES_OUTPUT_DIR && dir.exists(species_root)) {
    cat0("Removing old output directory: ", species_root)
    out_root_norm <- normalizePath(OUT_ROOT, winslash = "/", mustWork = TRUE)
    species_root_norm <- normalizePath(
      species_root,
      winslash = "/",
      mustWork = TRUE
    )
    expected_prefix <- paste0(out_root_norm, "/")
    if (!startsWith(species_root_norm, expected_prefix)) {
      stop(
        "Refusing to remove a SCENIC directory outside OUT_ROOT: ",
        species_root_norm
      )
    }
    unlink(species_root_norm, recursive = TRUE, force = TRUE)
  }

  out_dirs <- list(
    root = species_root,
    scenic = file.path(species_root, "SCENIC"),
    tables = file.path(species_root, "tables"),
    rds = file.path(species_root, "rds"),
    plots = file.path(species_root, "plots")
  )

  for (d in out_dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

  if (!file.exists(rds_file)) stop(label, ": input RDS not found: ", rds_file)

  check_db_files(db_dir, dbs, label)

  cat0("\n========== ", label, ": load object ==========")
  cat0("RDS: ", rds_file)
  obj <- readRDS(rds_file)
  if ("RNA" %in% Assays(obj)) DefaultAssay(obj) <- "RNA"

  cat0(label, " total cells in object: ", ncol(obj))
  cat0(label, " metadata columns:")
  print(colnames(obj@meta.data))

  obj <- make_status_fun(obj)

  prep <- prepare_expr_for_scenic(
    obj = obj,
    label = label,
    out_dirs = out_dirs,
    force_keep_genes = force_keep_genes,
    max_cells_per_group = MAX_CELLS_PER_GROUP,
    max_genes = MAX_GENES
  )

  rm(obj)
  gc()

  exprMat <- prep$exprMat
  cellInfo <- prep$cellInfo

  oldwd <- getwd()
  on.exit(setwd(oldwd), add = TRUE)
  setwd(out_dirs$scenic)

  dir.create("input", recursive = TRUE, showWarnings = FALSE)

  saveRDS(
    cellInfo,
    file = "input/cellInfo.Rds",
    compress = TRUE
  )

  write.csv(
    cellInfo,
    file = file.path(out_dirs$tables, paste0(label, "_cellInfo_for_SCENIC.csv")),
    row.names = FALSE,
    quote = FALSE
  )

  cat0("\n========== ", label, ": initializeScenic ==========")

  scenicOptions <- SCENIC::initializeScenic(
    org = org,
    dbDir = db_dir,
    dbs = dbs,
    datasetTitle = paste0(label, "_A4OL_formal_v9"),
    nCores = N_CORES,
    dbIndexCol = "motifs"
  )

  scenicOptions@settings$verbose <- TRUE
  scenicOptions@settings$nCores <- N_CORES
  scenicOptions@inputDatasetInfo$cellInfo <- "input/cellInfo.Rds"

  saveRDS(
    scenicOptions,
    file = file.path(out_dirs$rds, paste0(label, "_scenicOptions_init.rds")),
    compress = TRUE
  )

  cat0("\n========== ", label, ": runCorrelation ==========")
  SCENIC::runCorrelation(exprMat, scenicOptions)

  cat0("\n========== ", label, ": runGenie3 ==========")
  SCENIC::runGenie3(exprMat, scenicOptions)

  cat0("\n========== ", label, ": runSCENIC_1_coexNetwork2modules ==========")
  SCENIC::runSCENIC_1_coexNetwork2modules(scenicOptions)

  cat0("\n========== ", label, ": runSCENIC_2_createRegulons ==========")
  create_ok <- tryCatch({
    SCENIC::runSCENIC_2_createRegulons(
      scenicOptions,
      dbIndexCol = "motifs"
    )
    TRUE
  }, error = function(e) {
    cat0("[Warning] ", label, " createRegulons with dbIndexCol failed:")
    cat0(conditionMessage(e))
    cat0("[Retry] createRegulons without dbIndexCol, using patched importRankings.")
    SCENIC::runSCENIC_2_createRegulons(scenicOptions)
    TRUE
  })
  cat0(label, " createRegulons finished: ", create_ok)

  cat0("\n========== ", label, ": runSCENIC_3_scoreCells ==========")
  score_ok <- tryCatch({
    doMC::registerDoMC(cores = N_CORES)
    SCENIC::runSCENIC_3_scoreCells(
      scenicOptions,
      exprMat
    )
    TRUE
  }, error = function(e) {
    cat0("[Warning] ", label, " scoreCells stopped with message:")
    cat0(conditionMessage(e))
    cat0("If AUCell already finished and 3.4_regulonAUC.Rds exists, the script will continue extracting AUC.")
    FALSE
  })
  cat0(label, " scoreCells returned: ", score_ok)

  saveRDS(
    scenicOptions,
    file = file.path(out_dirs$rds, paste0(label, "_scenicOptions_finished_or_after_scoreCells.rds")),
    compress = TRUE
  )

  cat0("\n========== ", label, ": extract regulon AUC ==========")
  regulonAUC <- extract_auc_object(scenicOptions, out_dirs$scenic)
  auc_mat <- AUCell::getAUC(regulonAUC)

  cat0(label, " AUC matrix dim:")
  print(dim(auc_mat))

  write.csv(
    data.frame(regulon = rownames(auc_mat), TF = strip_regulon_tf(rownames(auc_mat))),
    file = file.path(out_dirs$tables, paste0(label, "_regulon_names.csv")),
    row.names = FALSE
  )

  write.csv(
    as.data.frame(as.matrix(auc_mat)),
    file = file.path(out_dirs$tables, paste0(label, "_regulon_AUC_matrix.csv")),
    quote = FALSE
  )

  saveRDS(
    regulonAUC,
    file = file.path(out_dirs$rds, paste0(label, "_regulonAUC.rds")),
    compress = TRUE
  )

  saveRDS(
    auc_mat,
    file = file.path(out_dirs$rds, paste0(label, "_auc_matrix.rds")),
    compress = TRUE
  )

  stats <- calc_regulon_stats(
    auc_mat = auc_mat,
    cellInfo = cellInfo,
    label = label,
    out_dirs = out_dirs
  )

  cat0("\n", label, " first 30 regulons:")
  print(head(rownames(auc_mat), 30))

  cat0("\n", label, " EGR1/Egr1 regulon stats:")
  print(stats %>% filter(TF_upper == "EGR1"))

  rm(exprMat, cellInfo, regulonAUC, auc_mat)
  gc()

  list(
    label = label,
    out_dirs = out_dirs,
    stats = stats,
    force_report = prep$force_report
  )
}

# ============================================================
# 3. Run human and mouse
# ============================================================
dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)

run_summary <- data.frame(
  parameter = c(
    "SEED",
    "N_CORES",
    "MAX_CELLS_PER_GROUP",
    "MAX_GENES",
    "MIN_CELLS_FACTOR",
    "MIN_COUNTS_FACTOR",
    "HUMAN_RDS",
    "MOUSE_RDS",
    "DB_DIR_HUMAN",
    "DB_DIR_MOUSE"
  ),
  value = c(
    SEED,
    N_CORES,
    MAX_CELLS_PER_GROUP,
    MAX_GENES,
    MIN_CELLS_FACTOR,
    MIN_COUNTS_FACTOR,
    HUMAN_RDS,
    MOUSE_RDS,
    DB_DIR_HUMAN,
    DB_DIR_MOUSE
  ),
  stringsAsFactors = FALSE
)

write.csv(run_summary, file = file.path(OUT_ROOT, "SCENIC_A4OL_FORMAL_run_parameters.csv"), row.names = FALSE)

human_res <- run_scenic_one_species(
  label = "human",
  org = "hgnc",
  rds_file = HUMAN_RDS,
  db_dir = DB_DIR_HUMAN,
  dbs = DBS_HUMAN,
  force_keep_genes = FORCE_KEEP_HUMAN,
  make_status_fun = make_human_a4_status
)

mouse_res <- run_scenic_one_species(
  label = "mouse",
  org = "mgi",
  rds_file = MOUSE_RDS,
  db_dir = DB_DIR_MOUSE,
  dbs = DBS_MOUSE,
  force_keep_genes = FORCE_KEEP_MOUSE,
  make_status_fun = make_mouse_a4_status
)

# ============================================================
# 4. Cross-species regulon comparison
# ============================================================
cat0("\n\n========== Cross-species regulon comparison ==========")

COMBINED_DIR <- file.path(OUT_ROOT, "combined_cross_species")
dir.create(COMBINED_DIR, recursive = TRUE, showWarnings = FALSE)

human_stats <- human_res$stats %>%
  mutate(species = "human") %>%
  select(
    species, regulon, TF, TF_upper,
    mean_other, mean_a4, median_other, median_a4,
    delta_mean_A4_minus_other, wilcox_p, wilcox_fdr
  )

mouse_stats <- mouse_res$stats %>%
  mutate(species = "mouse") %>%
  select(
    species, regulon, TF, TF_upper,
    mean_other, mean_a4, median_other, median_a4,
    delta_mean_A4_minus_other, wilcox_p, wilcox_fdr
  )

write.csv(
  bind_rows(human_stats, mouse_stats),
  file = file.path(COMBINED_DIR, "human_mouse_all_regulon_AUC_stats_long.csv"),
  row.names = FALSE
)

shared_tf <- inner_join(
  human_stats %>% rename_with(~ paste0("human_", .x), -TF_upper),
  mouse_stats %>% rename_with(~ paste0("mouse_", .x), -TF_upper),
  by = "TF_upper"
) %>%
  arrange(desc(human_delta_mean_A4_minus_other + mouse_delta_mean_A4_minus_other))

write.csv(
  shared_tf,
  file = file.path(COMBINED_DIR, "human_mouse_shared_regulon_TFs_wide.csv"),
  row.names = FALSE
)

shared_egr1 <- shared_tf %>% filter(TF_upper == "EGR1")
write.csv(
  shared_egr1,
  file = file.path(COMBINED_DIR, "human_mouse_EGR1_shared_regulon_summary.csv"),
  row.names = FALSE
)

force_combined <- bind_rows(
  human_res$force_report,
  mouse_res$force_report
)
write.csv(
  force_combined,
  file = file.path(COMBINED_DIR, "EGR1_Egr1_force_keep_gene_report_combined.csv"),
  row.names = FALSE
)

cat0("Shared regulon TF number: ", nrow(shared_tf))
cat0("EGR1 shared regulon summary:")
print(shared_egr1)

if (HAS_GGPLOT2 && nrow(shared_tf) > 0) {
  plot_df <- shared_tf %>%
    filter(is.finite(human_delta_mean_A4_minus_other), is.finite(mouse_delta_mean_A4_minus_other))

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = human_delta_mean_A4_minus_other, y = mouse_delta_mean_A4_minus_other)
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_point(size = 1.8, alpha = 0.75) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::labs(
      title = "Shared regulon TF activity change",
      x = "Human mean AUC difference: A4-OLs - Other OLs",
      y = "Mouse mean AUC difference: A4-OLs - Other OLs"
    )

  if ("EGR1" %in% plot_df$TF_upper) {
    p <- p + ggplot2::geom_point(
      data = plot_df %>% filter(TF_upper == "EGR1"),
      size = 3.2,
      shape = 21,
      stroke = 1.0
    )
  }

  ggplot2::ggsave(
    filename = file.path(COMBINED_DIR, "human_mouse_shared_regulon_delta_scatter.png"),
    plot = p,
    width = 5.5,
    height = 4.8,
    dpi = 300
  )

  ggplot2::ggsave(
    filename = file.path(COMBINED_DIR, "human_mouse_shared_regulon_delta_scatter.pdf"),
    plot = p,
    width = 5.5,
    height = 4.8
  )
}

cat0("\nDONE: formal human + mouse SCENIC A4-OLs analysis finished.")
cat0("Output root: ", OUT_ROOT)
cat0("Key outputs:")
cat0("  ", file.path(OUT_ROOT, "human/tables/human_all_regulon_AUC_A4_vs_Other_stats.csv"))
cat0("  ", file.path(OUT_ROOT, "mouse/tables/mouse_all_regulon_AUC_A4_vs_Other_stats.csv"))
cat0("  ", file.path(COMBINED_DIR, "human_mouse_shared_regulon_TFs_wide.csv"))
cat0("  ", file.path(COMBINED_DIR, "human_mouse_EGR1_shared_regulon_summary.csv"))
