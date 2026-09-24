#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# =============================================================================
# STEP5C_MOUSE_A4OL_LEAVE_ONE_MOUSE_OUT_ROBUSTNESS_FIXED_V4_ARROW_FIXED.R
#
# V4 修改：
#   1. 固定使用 oligo_annotated_cluster 作为原始四群注释来源
#      （Oligo1 / Oligo2 / Oligo3 / A4-OLs），不再错误识别 cell_type
#   2. 每轮重新聚类后，按原始四群的多数投票将新 cluster 合并为四群
#   3. UMAP 坐标轴改为 UMAP-1 / UMAP-2 箭头形式
#   4. 图中文字整体放大约两倍
#   5. 小鼠标签改为 A3_1~A3_3 和 A4_1~A4_5
#   6. 自动输出 8 轮三联图拼图，以及 8 轮四群重聚类 4×2 拼图
#
# V2 修复：
#   1. UCell 显式使用 name="_UCell"
#   2. 自动识别不同 UCell 版本生成的 score 列
#   3. 修复 RunHarmony 的 reduction 参数歧义
#   4. Harmony is applied by sample identity in every LOMO iteration
#
# 对每只 mouse 依次：
#   1) 完全删除该 mouse
#   2) 用剩余 mice 从 raw-count-derived data 重新标准化、降维、聚类
#   3) 只用剩余 mice 重新计算 A4-vs-Other markers 和 A4 signature
#   4) 用 UCell 重新评分
#   5) 自动识别 A4-like 新 cluster(s)
#   6) 重新计算 marker genes
#   7) 量化 signature recovery、marker overlap、UCell AUC 和 A4 enrichment
#
# 顺序运行，可断点续跑。每轮完成后写 DONE.txt。
# =============================================================================

options(stringsAsFactors = FALSE, timeout = 3600)
set.seed(20260727)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(UCell)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(showtext)
  library(sysfonts)
  library(patchwork)
})

source(file.path("R", "load_config.R"))

BASE_DIR <- A4OL_SERVER_ROOT
MOUSE_RDS <- file.path(
  A4OL_FIGURE_ROOT,
  "olig-c",
  "seurat_oligodendrocytes_subset_annotated_A4_OLs.rds"
)
OUT_ROOT <- file.path(
  A4OL_FIGURE_ROOT,
  "MOUSE_A4OL_validation_20260727",
  "03_leave_one_mouse_out_v4"
)
ITER_DIR <- file.path(OUT_ROOT, "iterations")
PLOT_DIR <- file.path(OUT_ROOT, "plots")
TABLE_DIR <- file.path(OUT_ROOT, "tables")
for (d in c(OUT_ROOT, ITER_DIR, PLOT_DIR, TABLE_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

SAMPLE_COL_OVERRIDE <- Sys.getenv("MOUSE_SAMPLE_COL", unset = "")

FORCE_RERUN <- tolower(
  Sys.getenv("LOMO_FORCE_RERUN", unset = "false")
) %in% c("true", "1", "yes")
# Harmony is applied using sample identity. APOE genotype is not supplied as a
# correction variable.
USE_HARMONY_IF_AVAILABLE <- TRUE
if (!requireNamespace("harmony", quietly = TRUE)) {
  stop("Package 'harmony' is required for the LOMO workflow.")
}
SAVE_LOMO_OBJECTS <- FALSE

TOP_SIGNATURE_N <- 200
TOP_RECLUSTER_MARKERS_N <- 200
MAX_CELLS_PER_IDENT_MARKERS <- 3000
N_VARIABLE_FEATURES <- 3000
N_PCS <- 30
CLUSTER_RESOLUTION <- 0.15
MIN_CLUSTER_CELLS <- 50
MAX_SELECTED_CLUSTERS <- 3
MIN_PCT <- 0.10
LOGFC_THRESHOLD <- 0.25

FONT_FILE <- A4OL_FONT_FILE
if (!file.exists(FONT_FILE)) stop("Arial font file not found: ", FONT_FILE)
sysfonts::font_add("Arial", regular = FONT_FILE)
showtext::showtext_auto()
showtext::showtext_opts(dpi = 300)
PLOT_FONT <- "Arial"

COL_A4 <- "#D19246"
COL_OTHER <- "#71A682"

# 原始四群配色：与正式图保持一致的柔绿色 + 橙色
FOUR_GROUP_COLORS <- c(
  "Oligo1" = "#B8DBB3",
  "Oligo2" = "#86BC79",
  "Oligo3" = "#71A682",
  "A4-OLs" = "#D19246"
)

FOUR_GROUP_LEVELS <- c(
  "Oligo1",
  "Oligo2",
  "Oligo3",
  "A4-OLs"
)

# 论文图使用的新动物名称
MOUSE_LABEL_MAP <- c(
  "A13" = "A3_1",
  "A14" = "A3_2",
  "A17" = "A3_3",
  "B17" = "A4_1",
  "B18" = "A4_2",
  "B19" = "A4_3",
  "B19right" = "A4_4",
  "B21" = "A4_5"
)

# 字体约为旧图的两倍
UMAP_TITLE_SIZE <- 26
UMAP_LEGEND_TITLE_SIZE <- 20
UMAP_LEGEND_TEXT_SIZE <- 18
UMAP_GROUP_LABEL_SIZE <- 9
UMAP_ARROW_TEXT_SIZE <- 5
SUMMARY_BASE_SIZE <- 24
SUMMARY_TITLE_SIZE <- 28
SUMMARY_STRIP_SIZE <- 22
SUMMARY_AXIS_TEXT_SIZE <- 19

cat0 <- function(...) cat(..., "\n", sep = "")

relabel_mouse <- function(x) {
  x <- as.character(x)
  out <- unname(MOUSE_LABEL_MAP[x])
  out[is.na(out)] <- x[is.na(out)]
  out
}

# 从输入对象建立原始四群。
# 优先使用 oligo_annotated_cluster；若不存在才使用原始 cluster 编号映射。
add_original_four_group <- function(obj) {
  meta_names <- as.character(names(obj@meta.data))

  if (any(meta_names == "oligo_annotated_cluster")) {
    x <- as.character(obj$oligo_annotated_cluster)

  } else if (any(meta_names == "oligo_clusters")) {
    cl <- as.character(obj$oligo_clusters)

    x <- dplyr::case_when(
      cl == "0" ~ "Oligo1",
      cl == "2" ~ "Oligo2",
      cl == "3" ~ "Oligo3",
      cl == "1" | cl == "4" ~ "A4-OLs",
      TRUE ~ NA_character_
    )

  } else {
    stop(
      "Cannot define the original four oligodendrocyte groups. ",
      "Required metadata: oligo_annotated_cluster or oligo_clusters."
    )
  }

  x <- trimws(x)

  # 兼容 Other/Other OLs：根据原始 cluster 补回 Oligo1/2/3。
  if (any(x %in% c("Other", "Other OLs"), na.rm = TRUE)) {
    if (!any(meta_names == "oligo_clusters")) {
      stop(
        "The annotation contains Other OLs but oligo_clusters is missing; ",
        "cannot separate Oligo1/Oligo2/Oligo3."
      )
    }

    cl <- as.character(obj$oligo_clusters)
    x[cl == "0"] <- "Oligo1"
    x[cl == "2"] <- "Oligo2"
    x[cl == "3"] <- "Oligo3"
    x[cl == "1" | cl == "4"] <- "A4-OLs"
  }

  unexpected <- setdiff(
    unique(x[!is.na(x)]),
    FOUR_GROUP_LEVELS
  )

  if (length(unexpected) > 0) {
    stop(
      "Unexpected values in the four-group annotation: ",
      paste(unexpected, collapse = ", ")
    )
  }

  obj$oligo_four_group <- factor(
    x,
    levels = FOUR_GROUP_LEVELS
  )

  if (any(is.na(obj$oligo_four_group))) {
    stop(
      "Some cells could not be assigned to the four oligodendrocyte groups."
    )
  }

  obj
}

# 每次重新聚类后，依据每个新 cluster 中原始四群数量最多的类别，
# 将新 cluster 合并为 Oligo1/Oligo2/Oligo3/A4-OLs。
map_reclustered_clusters_to_four_groups <- function(
  meta,
  cluster_col = "seurat_clusters"
) {
  map_table <- meta %>%
    tibble::rownames_to_column("cell") %>%
    transmute(
      cell = cell,
      new_cluster = as.character(.data[[cluster_col]]),
      original_four_group = as.character(oligo_four_group)
    ) %>%
    count(
      new_cluster,
      original_four_group,
      name = "n_cells"
    ) %>%
    group_by(new_cluster) %>%
    arrange(
      desc(n_cells),
      factor(
        original_four_group,
        levels = FOUR_GROUP_LEVELS
      ),
      .by_group = TRUE
    ) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    transmute(
      new_cluster = new_cluster,
      mapped_four_group = original_four_group
    )

  cluster_value <- as.character(meta[[cluster_col]])

  mapped <- map_table$mapped_four_group[
    match(
      cluster_value,
      map_table$new_cluster
    )
  ]

  list(
    group = factor(
      mapped,
      levels = FOUR_GROUP_LEVELS
    ),
    table = map_table
  )
}

# 将标准坐标轴替换为左下角箭头。
add_umap_arrow_axes <- function(
  p,
  obj,
  reduction = "umap"
) {
  emb <- as.data.frame(
    Seurat::Embeddings(
      obj,
      reduction = reduction
    )
  )

  xr <- range(emb[[1]], na.rm = TRUE)
  yr <- range(emb[[2]], na.rm = TRUE)

  dx <- diff(xr)
  dy <- diff(yr)

  # 将箭头整体移动到数据区域左下方，避免与 UMAP 点云重叠
  x0 <- xr[1] - 0.135 * dx
  y0 <- yr[1] - 0.135 * dy

  x1 <- x0 + 0.130 * dx
  y1 <- y0 + 0.130 * dy

  p +
    annotate(
      "segment",
      x = x0,
      xend = x1,
      y = y0,
      yend = y0,
      arrow = grid::arrow(
        length = grid::unit(0.22, "cm"),
        type = "closed"
      ),
      linewidth = 0.9,
      color = "black"
    ) +
    annotate(
      "segment",
      x = x0,
      xend = x0,
      y = y0,
      yend = y1,
      arrow = grid::arrow(
        length = grid::unit(0.22, "cm"),
        type = "closed"
      ),
      linewidth = 0.9,
      color = "black"
    ) +
    annotate(
      "text",
      x = (x0 + x1) / 2,
      y = y0 - 0.045 * dy,
      label = "UMAP-1",
      family = PLOT_FONT,
      size = UMAP_ARROW_TEXT_SIZE,
      color = "black"
    ) +
    annotate(
      "text",
      x = x0 - 0.045 * dx,
      y = (y0 + y1) / 2,
      label = "UMAP-2",
      angle = 90,
      family = PLOT_FONT,
      size = UMAP_ARROW_TEXT_SIZE,
      color = "black"
    ) +

    # 扩展左侧和下侧显示范围，仅为箭头预留空间；
    # UMAP 点的位置、颜色和其他图形设置均不改变
    coord_cartesian(
      xlim = c(
        xr[1] - 0.235 * dx,
        xr[2] + 0.015 * dx
      ),
      ylim = c(
        yr[1] - 0.235 * dy,
        yr[2] + 0.015 * dy
      ),
      clip = "off"
    ) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(
        t = 8,
        r = 8,
        b = 34,
        l = 34
      )
    )
}

large_umap_theme <- function() {
  theme_void(
    base_family = PLOT_FONT,
    base_size = 22
  ) +
    theme(
      text = element_text(
        family = PLOT_FONT,
        color = "black"
      ),
      plot.title = element_text(
        size = UMAP_TITLE_SIZE,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 8)
      ),
      legend.title = element_text(
        size = UMAP_LEGEND_TITLE_SIZE,
        face = "bold"
      ),
      legend.text = element_text(
        size = UMAP_LEGEND_TEXT_SIZE
      )
    )
}

pick_col <- function(meta, override = "", candidates, label, required = TRUE) {
  if (nzchar(override)) {
    if (!override %in% colnames(meta)) stop(label, " override not found: ", override)
    return(override)
  }
  hit <- candidates[candidates %in% colnames(meta)]
  if (length(hit) == 0) {
    if (required) stop("Cannot detect ", label)
    return(NA_character_)
  }
  hit[1]
}

add_a4_status <- function(obj) {
  meta <- obj@meta.data
  if ("oligo_annotated_cluster" %in% colnames(meta)) {
    obj$A4_status <- ifelse(
      as.character(obj$oligo_annotated_cluster) == "A4-OLs",
      "A4-OLs", "Other OLs"
    )
  } else if ("DEG_group" %in% colnames(meta)) {
    obj$A4_status <- ifelse(
      as.character(obj$DEG_group) == "A4-OLs",
      "A4-OLs", "Other OLs"
    )
  } else if ("oligo_clusters" %in% colnames(meta)) {
    obj$A4_status <- ifelse(
      as.character(obj$oligo_clusters) %in% c("1", "4"),
      "A4-OLs", "Other OLs"
    )
  } else {
    stop("Cannot define A4_status.")
  }
  obj$A4_status <- factor(obj$A4_status, levels = c("Other OLs", "A4-OLs"))
  obj
}

join_layers_safe <- function(obj) {
  has_join_layers <- exists(
    "JoinLayers",
    envir = asNamespace("SeuratObject"),
    mode = "function",
    inherits = FALSE
  )

  if (isTRUE(has_join_layers)) {
    obj <- tryCatch(
      SeuratObject::JoinLayers(
        object = obj,
        assay = "RNA"
      ),
      error = function(e) {
        cat0(
          "[JoinLayers warning] ",
          conditionMessage(e)
        )
        obj
      }
    )
  }

  obj
}


derive_signature <- function(obj, top_n, seed) {
  Idents(obj) <- "A4_status"
  mk <- FindMarkers(
    obj,
    ident.1 = "A4-OLs",
    ident.2 = "Other OLs",
    assay = "RNA",
    slot = "data",
    test.use = "wilcox",
    only.pos = TRUE,
    min.pct = MIN_PCT,
    logfc.threshold = LOGFC_THRESHOLD,
    max.cells.per.ident = MAX_CELLS_PER_IDENT_MARKERS,
    random.seed = seed,
    verbose = FALSE
  ) %>% rownames_to_column("gene")

  lfc_candidates <- c("avg_log2FC", "avg_logFC")
  lfc_col <- lfc_candidates[
    lfc_candidates %in% colnames(mk)
  ][1]

  if (is.na(lfc_col)) {
    stop("Cannot detect logFC column.")
  }

  mk <- mk %>%
    arrange(
      p_val_adj,
      desc(.data[[lfc_col]]),
      desc(pct.1 - pct.2)
    )

  list(
    markers = mk,
    genes = head(mk$gene, top_n)
  )
}

add_ucell_score_safe <- function(
  obj,
  signature_genes,
  signature_name = "LOMO_A4_signature"
) {
  signature_genes <- unique(
    intersect(
      as.character(signature_genes),
      rownames(obj)
    )
  )

  if (length(signature_genes) < 10) {
    stop(
      "Too few usable genes for UCell scoring: ",
      length(signature_genes)
    )
  }

  metadata_before <- colnames(obj@meta.data)

  feature_list <- list(signature_genes)
  names(feature_list) <- signature_name

  obj <- UCell::AddModuleScore_UCell(
    obj = obj,
    features = feature_list,
    assay = "RNA",
    name = "_UCell"
  )

  metadata_after <- colnames(obj@meta.data)
  expected_col <- paste0(
    signature_name,
    "_UCell"
  )

  # UCell 在当前环境中通常生成：
  # LOMO_A4_signature_UCell
  if (expected_col %in% metadata_after) {
    return(
      list(
        obj = obj,
        score_col = expected_col
      )
    )
  }

  # 兼容不同版本 UCell 的列名差异：
  # 优先从本次新增的数值列中识别。
  new_cols <- setdiff(
    metadata_after,
    metadata_before
  )

  numeric_new_cols <- new_cols[
    vapply(
      new_cols,
      function(x) {
        is.numeric(obj@meta.data[[x]])
      },
      logical(1)
    )
  ]

  score_candidates <- numeric_new_cols[
    grepl(
      paste0(
        signature_name,
        "|UCell|A4"
      ),
      numeric_new_cols,
      ignore.case = TRUE
    )
  ]

  # 若新增列没有命中，再在全部 metadata 中检索。
  if (length(score_candidates) == 0) {
    numeric_all_cols <- metadata_after[
      vapply(
        metadata_after,
        function(x) {
          is.numeric(obj@meta.data[[x]])
        },
        logical(1)
      )
    ]

    score_candidates <- numeric_all_cols[
      grepl(
        paste0(
          signature_name,
          "|LOMO.*UCell|A4.*UCell|UCell.*A4"
        ),
        numeric_all_cols,
        ignore.case = TRUE
      )
    ]
  }

  cat0(
    "[UCell] New metadata columns: ",
    ifelse(
      length(new_cols) == 0,
      "NONE",
      paste(new_cols, collapse = ", ")
    )
  )

  cat0(
    "[UCell] Candidate score columns: ",
    ifelse(
      length(score_candidates) == 0,
      "NONE",
      paste(score_candidates, collapse = ", ")
    )
  )

  if (length(score_candidates) == 0) {
    stop(
      "LOMO UCell score column missing. ",
      "Metadata columns after scoring: ",
      paste(metadata_after, collapse = ", ")
    )
  }

  list(
    obj = obj,
    score_col = score_candidates[1]
  )
}

manual_auc <- function(score, label) {
  ok <- is.finite(score) & !is.na(label)
  score <- score[ok]
  label <- as.logical(label[ok])
  n_pos <- sum(label)
  n_neg <- sum(!label)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  r <- rank(score, ties.method = "average")
  (sum(r[label]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

hypergeom_overlap_p <- function(overlap, set1, set2, universe) {
  if (any(c(overlap, set1, set2, universe) < 0) || universe <= 0) return(NA_real_)
  phyper(overlap - 1, set1, universe - set1, set2, lower.tail = FALSE)
}

cluster_enrichment_table <- function(meta, cluster_col, score_col) {
  global_frac <- mean(meta$A4_status == "A4-OLs")
  clusters <- sort(unique(as.character(meta[[cluster_col]])))

  out <- lapply(clusters, function(cl) {
    inside <- as.character(meta[[cluster_col]]) == cl
    a <- sum(inside & meta$A4_status == "A4-OLs")
    b <- sum(inside & meta$A4_status == "Other OLs")
    c <- sum(!inside & meta$A4_status == "A4-OLs")
    d <- sum(!inside & meta$A4_status == "Other OLs")

    ft <- tryCatch(
      fisher.test(
        matrix(c(a, b, c, d), nrow = 2),
        alternative = "greater"
      ),
      error = function(e) NULL
    )

    tibble(
      cluster = cl,
      n_cells = sum(inside),
      n_original_A4 = a,
      original_A4_fraction = ifelse(sum(inside) > 0, a / sum(inside), NA_real_),
      global_A4_fraction = global_frac,
      A4_enrichment_ratio = ifelse(
        global_frac > 0,
        original_A4_fraction / global_frac,
        NA_real_
      ),
      fisher_odds_ratio = ifelse(is.null(ft), NA_real_, unname(ft$estimate)),
      fisher_p = ifelse(is.null(ft), NA_real_, ft$p.value),
      median_A4_UCell = median(meta[[score_col]][inside], na.rm = TRUE),
      mean_A4_UCell = mean(meta[[score_col]][inside], na.rm = TRUE)
    )
  })

  bind_rows(out) %>%
    mutate(fisher_fdr = p.adjust(fisher_p, method = "BH")) %>%
    arrange(desc(median_A4_UCell), desc(A4_enrichment_ratio))
}

select_a4_clusters <- function(cluster_metrics) {
  eligible <- cluster_metrics %>% filter(n_cells >= MIN_CLUSTER_CELLS)
  if (nrow(eligible) == 0) eligible <- cluster_metrics

  score_q75 <- quantile(
    eligible$median_A4_UCell,
    0.75,
    na.rm = TRUE
  )

  selected <- eligible %>%
    filter(
      median_A4_UCell >= score_q75,
      A4_enrichment_ratio > 1,
      fisher_p < 0.05
    ) %>%
    arrange(desc(median_A4_UCell), desc(A4_enrichment_ratio))

  if (nrow(selected) > MAX_SELECTED_CLUSTERS) {
    selected <- selected[seq_len(MAX_SELECTED_CLUSTERS), , drop = FALSE]
  }

  if (nrow(selected) == 0) {
    selected <- eligible %>%
      arrange(desc(median_A4_UCell), desc(A4_enrichment_ratio))
    selected <- selected[1, , drop = FALSE]
  }

  selected$cluster
}

run_reclustering <- function(obj, batch_col, seed) {
  set.seed(seed)
  DefaultAssay(obj) <- "RNA"

  obj <- NormalizeData(
    obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )
  obj <- FindVariableFeatures(
    obj,
    selection.method = "vst",
    nfeatures = N_VARIABLE_FEATURES,
    verbose = FALSE
  )
  obj <- ScaleData(
    obj,
    features = VariableFeatures(obj),
    verbose = FALSE
  )
  obj <- RunPCA(
    obj,
    features = VariableFeatures(obj),
    npcs = N_PCS,
    seed.use = seed,
    verbose = FALSE
  )

  reduction_use <- "pca"

  if (
    USE_HARMONY_IF_AVAILABLE &&
    requireNamespace("harmony", quietly = TRUE) &&
    !is.na(batch_col) &&
    batch_col %in% colnames(obj@meta.data) &&
    n_distinct(obj@meta.data[[batch_col]]) > 1
  ) {
    obj_h <- tryCatch(
      harmony::RunHarmony(
        object = obj,
        group.by.vars = batch_col,
        reduction.use = "pca",
        dims.use = seq_len(
          min(
            N_PCS,
            ncol(
              Embeddings(
                obj,
                reduction = "pca"
              )
            )
          )
        ),
        assay.use = "RNA",
        project.dim = FALSE,
        verbose = FALSE
      ),
      error = function(e) {
        cat0("[Harmony warning] ", conditionMessage(e))
        NULL
      }
    )

    if (!is.null(obj_h) && "harmony" %in% Reductions(obj_h)) {
      obj <- obj_h
      reduction_use <- "harmony"
    }
  }

  nd <- min(
    N_PCS,
    ncol(Embeddings(obj, reduction = reduction_use))
  )

  obj <- FindNeighbors(
    obj,
    reduction = reduction_use,
    dims = seq_len(nd),
    verbose = FALSE
  )
  obj <- FindClusters(
    obj,
    resolution = CLUSTER_RESOLUTION,
    random.seed = seed,
    verbose = FALSE
  )
  obj <- RunUMAP(
    obj,
    reduction = reduction_use,
    dims = seq_len(nd),
    seed.use = seed,
    verbose = FALSE
  )

  list(obj = obj, reduction = reduction_use)
}

if (!file.exists(MOUSE_RDS)) stop("Mouse RDS not found: ", MOUSE_RDS)
obj0 <- readRDS(MOUSE_RDS)

assay_names <- as.character(
  names(obj0@assays)
)

if (!any(assay_names == "RNA")) {
  stop(
    "RNA assay missing. Available assays: ",
    paste(assay_names, collapse = ", ")
  )
}

DefaultAssay(obj0) <- "RNA"
obj0 <- join_layers_safe(obj0)
obj0 <- add_a4_status(obj0)
obj0 <- add_original_four_group(obj0)

meta0 <- obj0@meta.data
sample_col <- pick_col(
  meta0,
  SAMPLE_COL_OVERRIDE,
  c("sample", "original_sample", "orig.ident", "sample_id", "donor"),
  "mouse/sample column"
)
# Use biological sample identity for Harmony integration in each LOMO run.
# APOE genotype is deliberately not used as a correction variable.
batch_col <- sample_col

obj0$mouse_id <- as.character(obj0@meta.data[[sample_col]])
obj0$mouse_plot_label <- relabel_mouse(obj0$mouse_id)

mice <- sort(unique(
  obj0$mouse_id[!is.na(obj0$mouse_id) & obj0$mouse_id != ""]
))

mouse_plot_levels <- relabel_mouse(mice)

cat0("Mouse ID column: ", sample_col)
cat0("Batch column: ", ifelse(is.na(batch_col), "NONE", batch_col))
cat0(
  "Harmony enabled: ",
  USE_HARMONY_IF_AVAILABLE
)
cat0("Number of mice: ", length(mice))
cat0("\nOriginal four-group distribution:")
print(
  table(
    relabel_mouse(obj0$mouse_id),
    obj0$oligo_four_group
  )
)

cat0("\nBinary A4 status:")
print(
  table(
    relabel_mouse(obj0$mouse_id),
    obj0$A4_status
  )
)

obj0 <- NormalizeData(obj0, verbose = FALSE)
baseline <- derive_signature(obj0, TOP_SIGNATURE_N, 20260727)
baseline_genes <- baseline$genes

write_csv(
  baseline$markers,
  file.path(TABLE_DIR, "baseline_A4_vs_Other_positive_markers.csv")
)
write_csv(
  tibble(rank = seq_along(baseline_genes), gene = baseline_genes),
  file.path(TABLE_DIR, "baseline_A4_signature_top200.csv")
)

summary_rows <- list()
marker_rows <- list()

# 保存每轮三联图及四群重聚类图，用于最后自动拼图
iteration_triptych_list <- list()
iteration_four_group_list <- list()

for (i in seq_along(mice)) {
  held_out <- mice[i]
  held_out_plot <- relabel_mouse(held_out)
  safe_name <- gsub("[^A-Za-z0-9_.-]", "_", held_out_plot)
  iter_out <- file.path(ITER_DIR, paste0("leave_out_", safe_name))
  dir.create(iter_out, recursive = TRUE, showWarnings = FALSE)

  done_file <- file.path(iter_out, "DONE.txt")
  summary_file <- file.path(iter_out, "iteration_summary.csv")

  if (file.exists(done_file) &&
      file.exists(summary_file) &&
      !FORCE_RERUN) {
    cat0("[Reuse] ", held_out)
    summary_rows[[length(summary_rows) + 1]] <-
      read.csv(summary_file, check.names = FALSE)

    marker_file <- file.path(iter_out, "reclustered_A4_markers_top200.csv")
    if (file.exists(marker_file)) {
      x <- read.csv(marker_file, check.names = FALSE)
      x$held_out_mouse <- held_out
      marker_rows[[length(marker_rows) + 1]] <- x
    }
    next
  }

  cat0("\n============================================================")
  cat0("LOMO ", i, "/", length(mice), ": leave out ", held_out_plot, " [", held_out, "]")
  cat0("============================================================")

  cells_train <- colnames(obj0)[obj0$mouse_id != held_out]
  train <- subset(obj0, cells = cells_train)
  train <- join_layers_safe(train)
  train <- NormalizeData(train, verbose = FALSE)

  # signature 完全由剩余 mice 生成
  sig_res <- derive_signature(
    train,
    TOP_SIGNATURE_N,
    20260727 + i
  )
  train_signature <- sig_res$genes

  write_csv(
    sig_res$markers,
    file.path(iter_out, "training_A4_vs_Other_positive_markers.csv")
  )
  write_csv(
    tibble(rank = seq_along(train_signature), gene = train_signature),
    file.path(iter_out, "training_A4_signature_top200.csv")
  )

  # 从头重新 clustering
  recl <- run_reclustering(
    train,
    batch_col,
    20260727 + i
  )
  train <- recl$obj

  # UCell：显式指定 name="_UCell"，并兼容不同 UCell 版本的列名
  ucell_res <- add_ucell_score_safe(
    obj = train,
    signature_genes = train_signature,
    signature_name = "LOMO_A4_signature"
  )

  train <- ucell_res$obj
  score_col <- ucell_res$score_col

  cat0(
    "[OK] LOMO UCell score column: ",
    score_col
  )

  cat0(
    "[OK] LOMO UCell score summary:"
  )
  print(
    summary(
      train@meta.data[[score_col]]
    )
  )

  cluster_col <- "seurat_clusters"
  cm <- cluster_enrichment_table(
    train@meta.data,
    cluster_col,
    score_col
  )
  selected_clusters <- select_a4_clusters(cm)
  cm$selected_as_A4 <- cm$cluster %in% selected_clusters
  write_csv(
    cm,
    file.path(iter_out, "reclustered_cluster_A4_metrics.csv")
  )

  train$LOMO_A4_call <- ifelse(
    as.character(train@meta.data[[cluster_col]]) %in% selected_clusters,
    "LOMO A4-like cluster",
    "Other clusters"
  )
  train$LOMO_A4_call <- factor(
    train$LOMO_A4_call,
    levels = c("Other clusters", "LOMO A4-like cluster")
  )

  # 新聚类按照原始四群多数投票合并为四个群
  four_map <- map_reclustered_clusters_to_four_groups(
    meta = train@meta.data,
    cluster_col = cluster_col
  )

  train$LOMO_four_group <- four_map$group

  write_csv(
    four_map$table,
    file.path(
      iter_out,
      "reclustered_cluster_to_four_group_mapping.csv"
    )
  )

  cat0("Reclustered four-group distribution:")
  print(
    table(
      train$LOMO_four_group,
      useNA = "ifany"
    )
  )

  Idents(train) <- "LOMO_A4_call"
  recluster_markers <- FindMarkers(
    train,
    ident.1 = "LOMO A4-like cluster",
    ident.2 = "Other clusters",
    assay = "RNA",
    slot = "data",
    test.use = "wilcox",
    only.pos = TRUE,
    min.pct = MIN_PCT,
    logfc.threshold = LOGFC_THRESHOLD,
    max.cells.per.ident = MAX_CELLS_PER_IDENT_MARKERS,
    random.seed = 20260727 + i,
    verbose = FALSE
  ) %>% rownames_to_column("gene")

  lfc_col <- c("avg_log2FC", "avg_logFC")[
    c("avg_log2FC", "avg_logFC") %in% colnames(recluster_markers)
  ][1]
  if (is.na(lfc_col)) stop("Cannot detect recluster marker logFC column.")

  recluster_markers <- recluster_markers %>%
    arrange(p_val_adj, desc(.data[[lfc_col]]), desc(pct.1 - pct.2))

  n_keep <- min(TOP_RECLUSTER_MARKERS_N, nrow(recluster_markers))
  top_recluster <- recluster_markers[
    seq_len(n_keep),
    ,
    drop = FALSE
  ]

  write_csv(
    recluster_markers,
    file.path(iter_out, "reclustered_A4_markers_all.csv")
  )
  write_csv(
    top_recluster,
    file.path(iter_out, "reclustered_A4_markers_top200.csv")
  )

  marker_rows[[length(marker_rows) + 1]] <- top_recluster %>%
    mutate(held_out_mouse = held_out_plot, held_out_mouse_raw = held_out)

  overlap_genes <- intersect(top_recluster$gene, baseline_genes)
  universe_n <- length(rownames(train))
  overlap_p <- hypergeom_overlap_p(
    length(overlap_genes),
    length(baseline_genes),
    nrow(top_recluster),
    universe_n
  )

  meta <- train@meta.data
  call <- meta$LOMO_A4_call == "LOMO A4-like cluster"
  score <- meta[[score_col]]

  score_delta <- median(score[call], na.rm = TRUE) -
    median(score[!call], na.rm = TRUE)
  score_p <- tryCatch(
    wilcox.test(score[call], score[!call])$p.value,
    error = function(e) NA_real_
  )
  auc <- manual_auc(score, call)

  selected_metric <- cm %>% filter(selected_as_A4)

  iter_summary <- tibble(
    held_out_mouse = held_out_plot,
    held_out_mouse_raw = held_out,
    n_training_mice = n_distinct(train$mouse_id),
    n_training_cells = ncol(train),
    n_training_original_A4_cells = sum(train$A4_status == "A4-OLs"),
    n_reclustered_clusters = n_distinct(train$seurat_clusters),
    selected_clusters = paste(selected_clusters, collapse = ";"),
    n_selected_clusters = length(selected_clusters),
    selected_cluster_cells = sum(call),
    selected_cluster_original_A4_fraction =
      mean(train$A4_status[call] == "A4-OLs"),
    selected_cluster_median_UCell = median(score[call], na.rm = TRUE),
    other_cluster_median_UCell = median(score[!call], na.rm = TRUE),
    UCell_median_delta = score_delta,
    UCell_wilcox_p = score_p,
    UCell_ROC_AUC = auc,
    marker_overlap_n = length(overlap_genes),
    marker_recovery_fraction =
      length(overlap_genes) / length(baseline_genes),
    marker_jaccard =
      length(overlap_genes) /
      length(union(top_recluster$gene, baseline_genes)),
    marker_overlap_hypergeom_p = overlap_p,
    max_A4_enrichment_ratio =
      max(selected_metric$A4_enrichment_ratio, na.rm = TRUE),
    min_cluster_fisher_p =
      min(selected_metric$fisher_p, na.rm = TRUE),
    reduction_used = recl$reduction
  )

  write_csv(iter_summary, summary_file)
  writeLines(
    overlap_genes,
    file.path(iter_out, "baseline_signature_overlap_genes.txt")
  )

  # ---------------------------------------------------------------------------
  # 四群重聚类图：新 cluster 按原始四群多数投票合并
  # ---------------------------------------------------------------------------
  p_cluster <- DimPlot(
    train,
    reduction = "umap",
    group.by = "LOMO_four_group",
    label = TRUE,
    repel = TRUE,
    label.size = UMAP_GROUP_LABEL_SIZE,
    pt.size = 0.55,
    cols = FOUR_GROUP_COLORS,
    raster = FALSE
  ) +
    ggtitle(
      paste0(
        "Leave out ",
        held_out_plot,
        ": four-group reclustering"
      )
    ) +
    large_umap_theme() +
    theme(
      legend.position = "none"
    )

  p_cluster <- add_umap_arrow_axes(
    p = p_cluster,
    obj = train,
    reduction = "umap"
  )

  # ---------------------------------------------------------------------------
  # Training-only UCell score
  # ---------------------------------------------------------------------------
  p_score <- FeaturePlot(
    train,
    features = score_col,
    reduction = "umap",
    order = TRUE,
    min.cutoff = "q05",
    max.cutoff = "q95",
    pt.size = 0.55,
    raster = FALSE
  ) +
    scale_color_gradientn(
      colours = c(
        "#F0EEF6",
        "#B9B2D8",
        "#6958B5",
        "#271C78"
      ),
      name = "UCell\nscore"
    ) +
    ggtitle(
      "Training-only A4 UCell score"
    ) +
    large_umap_theme() +
    theme(
      legend.position = "right",
      legend.key.height = grid::unit(1.4, "cm"),
      legend.key.width = grid::unit(0.55, "cm")
    )

  p_score <- add_umap_arrow_axes(
    p = p_score,
    obj = train,
    reduction = "umap"
  )

  # ---------------------------------------------------------------------------
  # 自动识别的 A4-like cluster
  # ---------------------------------------------------------------------------
  p_call <- DimPlot(
    train,
    reduction = "umap",
    group.by = "LOMO_A4_call",
    pt.size = 0.55,
    raster = FALSE,
    cols = c(
      "Other clusters" = "grey82",
      "LOMO A4-like cluster" = COL_A4
    )
  ) +
    ggtitle(
      "Automatically rediscovered A4-like cluster(s)"
    ) +
    large_umap_theme() +
    theme(
      legend.position = "right",
      legend.title = element_blank()
    )

  p_call <- add_umap_arrow_axes(
    p = p_call,
    obj = train,
    reduction = "umap"
  )

  p_combined <- p_cluster + p_score + p_call +
    plot_layout(
      ncol = 3,
      widths = c(1, 1.08, 1.08)
    ) +
    plot_annotation(
      title = paste0(
        "Leave-one-mouse-out: ",
        held_out_plot
      ),
      theme = theme(
        text = element_text(
          family = PLOT_FONT
        ),
        plot.title = element_text(
          face = "bold",
          hjust = 0.5,
          size = 30,
          margin = margin(b = 12)
        )
      )
    )

  ggsave(
    file.path(
      iter_out,
      "LOMO_reclustering_A4_signature.png"
    ),
    p_combined,
    width = 24,
    height = 7.8,
    dpi = 300,
    bg = "white",
    limitsize = FALSE
  )

  ggsave(
    file.path(
      iter_out,
      "LOMO_reclustering_A4_signature.pdf"
    ),
    p_combined,
    width = 24,
    height = 7.8,
    bg = "white",
    limitsize = FALSE
  )

  iteration_triptych_list[[held_out_plot]] <- p_combined
  iteration_four_group_list[[held_out_plot]] <- p_cluster

  if (SAVE_LOMO_OBJECTS) {
    saveRDS(
      train,
      file.path(iter_out, "LOMO_reclustered_object.rds"),
      compress = TRUE
    )
  }

  writeLines(as.character(Sys.time()), done_file)
  summary_rows[[length(summary_rows) + 1]] <- iter_summary

  rm(train, recl, sig_res, recluster_markers)
  gc()
}

# =============================================================================
# 拼接全部 8 轮 LOMO 图
# =============================================================================
if (length(iteration_triptych_list) > 0) {
  combined_triptych <- patchwork::wrap_plots(
    iteration_triptych_list,
    ncol = 1
  ) +
    plot_annotation(
      title =
        "Leave-one-mouse-out A4-OL reclustering and signature recovery",
      theme = theme(
        text = element_text(
          family = PLOT_FONT
        ),
        plot.title = element_text(
          face = "bold",
          hjust = 0.5,
          size = 32,
          margin = margin(b = 16)
        )
      )
    )

  ggsave(
    file.path(
      PLOT_DIR,
      "LOMO_reclustering_A4_signature_all_8_combined.png"
    ),
    combined_triptych,
    width = 24,
    height = 7.9 * length(iteration_triptych_list),
    dpi = 300,
    bg = "white",
    limitsize = FALSE
  )

  ggsave(
    file.path(
      PLOT_DIR,
      "LOMO_reclustering_A4_signature_all_8_combined.pdf"
    ),
    combined_triptych,
    width = 24,
    height = 7.9 * length(iteration_triptych_list),
    bg = "white",
    limitsize = FALSE
  )
}

# 只拼接每轮左侧四群重聚类图，4×2 排版，适合作为文章 panel
if (length(iteration_four_group_list) > 0) {
  combined_four_group <- patchwork::wrap_plots(
    iteration_four_group_list,
    ncol = 4
  ) +
    plot_annotation(
      title =
        "Four-group reclustering across leave-one-mouse-out iterations",
      theme = theme(
        text = element_text(
          family = PLOT_FONT
        ),
        plot.title = element_text(
          face = "bold",
          hjust = 0.5,
          size = 32,
          margin = margin(b = 14)
        )
      )
    )

  ggsave(
    file.path(
      PLOT_DIR,
      "LOMO_four_group_reclustering_all_8_4x2.png"
    ),
    combined_four_group,
    width = 24,
    height = 13,
    dpi = 300,
    bg = "white",
    limitsize = FALSE
  )

  ggsave(
    file.path(
      PLOT_DIR,
      "LOMO_four_group_reclustering_all_8_4x2.pdf"
    ),
    combined_four_group,
    width = 24,
    height = 13,
    bg = "white",
    limitsize = FALSE
  )
}

summary_df <- bind_rows(summary_rows) %>%
  mutate(
    UCell_wilcox_fdr = p.adjust(UCell_wilcox_p, method = "BH"),
    marker_overlap_hypergeom_fdr =
      p.adjust(marker_overlap_hypergeom_p, method = "BH")
  )
write_csv(
  summary_df,
  file.path(TABLE_DIR, "LOMO_all_iterations_summary.csv")
)

markers_all <- bind_rows(marker_rows)
if (nrow(markers_all) > 0) {
  recurrence <- markers_all %>%
    distinct(held_out_mouse, gene) %>%
    count(gene, name = "n_LOMO_iterations") %>%
    mutate(
      total_iterations = n_distinct(summary_df$held_out_mouse),
      recurrence_fraction = n_LOMO_iterations / total_iterations,
      is_baseline_signature = gene %in% baseline_genes
    ) %>%
    arrange(
      desc(n_LOMO_iterations),
      desc(is_baseline_signature),
      gene
    )

  write_csv(
    recurrence,
    file.path(TABLE_DIR, "LOMO_marker_recurrence.csv")
  )

  top_rec <- recurrence[
    seq_len(min(30, nrow(recurrence))),
    ,
    drop = FALSE
  ]
  top_rec$gene <- factor(top_rec$gene, levels = rev(top_rec$gene))

  p_rec <- ggplot(
    top_rec,
    aes(n_LOMO_iterations, gene, color = is_baseline_signature)
  ) +
    geom_segment(
      aes(
        x = 0,
        xend = n_LOMO_iterations,
        y = gene,
        yend = gene
      ),
      color = "grey75",
      linewidth = 0.5
    ) +
    geom_point(size = 3.2) +
    scale_color_manual(values = c(
      "TRUE" = "#D19246",
      "FALSE" = "grey45"
    )) +
    labs(
      title = "A4 marker recurrence across leave-one-mouse-out analyses",
      x = "Number of LOMO iterations",
      y = NULL,
      color = "Baseline A4\nsignature"
    ) +
    theme_classic(base_family = PLOT_FONT, base_size = 22) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5)
    )

  ggsave(
    file.path(PLOT_DIR, "LOMO_A4_marker_recurrence_top30.png"),
    p_rec,
    width = 7.5,
    height = 7.5,
    dpi = 300,
    bg = "white"
  )
  ggsave(
    file.path(PLOT_DIR, "LOMO_A4_marker_recurrence_top30.pdf"),
    p_rec,
    width = 7.5,
    height = 7.5,
    bg = "white"
  )
}

summary_long <- summary_df %>%
  select(
    held_out_mouse,
    marker_recovery_fraction,
    UCell_ROC_AUC,
    UCell_median_delta,
    selected_cluster_original_A4_fraction
  ) %>%
  pivot_longer(
    -held_out_mouse,
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(
      metric,
      marker_recovery_fraction = "Baseline marker recovery",
      UCell_ROC_AUC = "UCell ROC AUC",
      UCell_median_delta = "UCell median delta",
      selected_cluster_original_A4_fraction =
        "Original A4 fraction in rediscovered cluster"
    )
  )

p_summary <- ggplot(
  summary_long,
  aes(held_out_mouse, value)
) +
  geom_col(fill = "#D19246", alpha = 0.85, width = 0.75) +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  labs(
    title = "Leave-one-mouse-out robustness of the A4-OL signature",
    x = "Held-out mouse",
    y = NULL
  ) +
  theme_classic(
    base_family = PLOT_FONT,
    base_size = SUMMARY_BASE_SIZE
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = SUMMARY_TITLE_SIZE
    ),
    axis.title.x = element_text(
      size = SUMMARY_BASE_SIZE
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = SUMMARY_AXIS_TEXT_SIZE,
      color = "black"
    ),
    axis.text.y = element_text(
      size = SUMMARY_AXIS_TEXT_SIZE,
      color = "black"
    ),
    strip.text = element_text(
      face = "bold",
      size = SUMMARY_STRIP_SIZE
    ),
    strip.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.8
    )
  )

ggsave(
  file.path(PLOT_DIR, "LOMO_A4_signature_robustness_summary.png"),
  p_summary,
  width = 16,
  height = 11,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(PLOT_DIR, "LOMO_A4_signature_robustness_summary.pdf"),
  p_summary,
  width = 16,
  height = 11,
  bg = "white"
)

cat0("\nDONE")
cat0("Output: ", OUT_ROOT)
