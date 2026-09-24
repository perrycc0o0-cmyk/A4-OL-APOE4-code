#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# =============================================================================
# STEP5A_MOUSE_A4OL_UCELL_APOE3_VS_APOE4_WITH_SIGNIFICANCE.R
#
# 目的：
#   1. 在小鼠 oligodendrocyte 对象中读取或构建固定 A4-OL signature
#   2. 使用 UCell 为每个细胞计算 A4-OL signature score
#   3. 在 A4-OL 亚群内比较 APOE3 vs APOE4
#   4. 输出 cell-level 和 mouse-level 统计
#   5. 在两张图上添加显著性括号和星号
#
# 统计原则：
#   - cell-level Wilcoxon 仅作描述，因为同一只 mouse 内的细胞并非独立重复
#   - mouse-level Wilcoxon 为主要统计结果，每个点代表一只独立 mouse
#   - B19 和 B19right 按两个独立样本分别处理
#   - 排除 A4-OL 细胞数少于 10 的 mouse（本数据中为 A17、B17）
# =============================================================================

options(stringsAsFactors = FALSE)
options(timeout = 1200)
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
})

source(file.path("R", "load_config.R"))

# =============================================================================
# 0. 路径和参数
# =============================================================================
BASE_DIR <- A4OL_SERVER_ROOT

MOUSE_RDS <- file.path(
  A4OL_FIGURE_ROOT,
  "olig-c",
  "seurat_oligodendrocytes_subset_annotated_A4_OLs.rds"
)

OUT_ROOT <- file.path(
  A4OL_FIGURE_ROOT,
  "MOUSE_A4OL_validation_20260727",
  "01_UCell_APOE3_vs_APOE4"
)

PLOT_DIR <- file.path(OUT_ROOT, "plots")
TABLE_DIR <- file.path(OUT_ROOT, "tables")
RDS_DIR <- file.path(OUT_ROOT, "rds")

for (d in c(OUT_ROOT, PLOT_DIR, TABLE_DIR, RDS_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# 若已有固定 A4 signature，推荐运行前指定：
# export A4_SIGNATURE_FILE=/path/to/mouse_A4_signature.tsv
#
# 支持 csv/tsv/txt，脚本会自动识别 gene 列。
SIGNATURE_FILE <- Sys.getenv("A4_SIGNATURE_FILE", unset = "")

# 自动构建 signature 时使用 A4-OL vs Other OLs 的前多少个上调 marker
TOP_SIGNATURE_N <- 200L

# FindMarkers 参数
MIN_PCT <- 0.10
LOGFC_THRESHOLD <- 0.25
MAX_CELLS_PER_IDENT_MARKERS <- 3000L

# 每只 mouse 至少需要 10 个 A4-OL 才纳入后续比较
# 因此会自动排除 A17 和 B17（各只有 1 个 A4-OL）
MIN_A4_CELLS_PER_MOUSE <- 10L

# metadata 列覆盖。通常不需要设置。
# export MOUSE_SAMPLE_COL=sample
# export MOUSE_GENOTYPE_COL=genotype
SAMPLE_COL_OVERRIDE <- Sys.getenv("MOUSE_SAMPLE_COL", unset = "")
GENOTYPE_COL_OVERRIDE <- Sys.getenv("MOUSE_GENOTYPE_COL", unset = "")

# Arial 字体
FONT_FILE <- A4OL_FONT_FILE

if (!file.exists(FONT_FILE)) {
  stop("Arial font file not found: ", FONT_FILE)
}

sysfonts::font_add("Arial", regular = FONT_FILE)
showtext::showtext_auto()
showtext::showtext_opts(dpi = 300)

PLOT_FONT <- "Arial"

COLORS <- c(
  "APOE3" = "#4DBBD5",
  "APOE4" = "#E64B35"
)

# 图形参数
CELL_PLOT_WIDTH <- 9
CELL_PLOT_HEIGHT <- 6.0
MOUSE_PLOT_WIDTH <- 6.4
MOUSE_PLOT_HEIGHT <- 6.0

TITLE_SIZE <- 18
SUBTITLE_SIZE <- 13
AXIS_TITLE_SIZE <- 16
AXIS_TEXT_SIZE <- 14
SIG_TEXT_SIZE <- 6
SIG_LINE_WIDTH <- 0.75

# =============================================================================
# 1. 通用函数
# =============================================================================
cat0 <- function(...) {
  cat(..., "\n", sep = "")
}

pick_col <- function(
  meta,
  override = "",
  candidates,
  label
) {
  if (nzchar(override)) {
    if (!override %in% colnames(meta)) {
      stop(label, " override column not found: ", override)
    }
    return(override)
  }

  hit <- candidates[candidates %in% colnames(meta)]

  if (length(hit) == 0) {
    stop(
      "Cannot detect ",
      label,
      ". Candidate columns: ",
      paste(candidates, collapse = ", ")
    )
  }

  hit[1]
}

recode_genotype <- function(x) {
  z <- toupper(
    gsub(
      "[^A-Z0-9]",
      "",
      as.character(x)
    )
  )

  dplyr::case_when(
    z %in% c(
      "APOE3",
      "E3",
      "3",
      "E3E3",
      "APOE33",
      "33"
    ) ~ "APOE3",

    z %in% c(
      "APOE4",
      "E4",
      "4",
      "E4E4",
      "APOE44",
      "44",
      "E3E4",
      "APOE34",
      "34"
    ) ~ "APOE4",

    grepl("APOE4|E4", z) ~ "APOE4",
    grepl("APOE3|E3", z) ~ "APOE3",
    TRUE ~ NA_character_
  )
}

add_a4_status <- function(obj) {
  meta <- obj@meta.data

  if ("oligo_annotated_cluster" %in% colnames(meta)) {
    obj$A4_status <- ifelse(
      as.character(obj$oligo_annotated_cluster) == "A4-OLs",
      "A4-OLs",
      "Other OLs"
    )

  } else if ("DEG_group" %in% colnames(meta)) {
    obj$A4_status <- ifelse(
      as.character(obj$DEG_group) == "A4-OLs",
      "A4-OLs",
      "Other OLs"
    )

  } else if ("oligo_clusters" %in% colnames(meta)) {
    obj$A4_status <- ifelse(
      as.character(obj$oligo_clusters) %in% c("1", "4"),
      "A4-OLs",
      "Other OLs"
    )

  } else {
    stop(
      "Cannot define A4_status. Missing ",
      "oligo_annotated_cluster / DEG_group / oligo_clusters."
    )
  }

  obj$A4_status <- factor(
    obj$A4_status,
    levels = c("Other OLs", "A4-OLs")
  )

  obj
}

join_layers_safe <- function(obj) {
  if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
    obj <- tryCatch(
      SeuratObject::JoinLayers(
        obj,
        assay = "RNA"
      ),
      error = function(e) {
        cat0("[Warning] JoinLayers failed: ", conditionMessage(e))
        obj
      }
    )
  }

  obj
}

read_signature <- function(
  path,
  available_genes
) {
  ext <- tolower(
    tools::file_ext(path)
  )

  if (ext == "csv") {
    x <- read.csv(
      path,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

  } else if (ext %in% c("tsv", "txt")) {
    x <- read.delim(
      path,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

  } else {
    stop(
      "Unsupported signature file extension: ",
      ext
    )
  }

  candidates <- c(
    "gene",
    "Gene",
    "gene_symbol",
    "symbol",
    "features",
    "target_gene"
  )

  gene_col <- candidates[
    candidates %in% colnames(x)
  ][1]

  if (is.na(gene_col)) {
    gene_col <- colnames(x)[1]
  }

  genes <- unique(
    as.character(
      x[[gene_col]]
    )
  )

  genes <- genes[
    !is.na(genes) &
      genes != ""
  ]

  intersect(
    genes,
    available_genes
  )
}

derive_a4_signature <- function(
  obj,
  top_n = 200L
) {
  Idents(obj) <- "A4_status"

  markers <- FindMarkers(
    object = obj,
    ident.1 = "A4-OLs",
    ident.2 = "Other OLs",
    assay = "RNA",
    slot = "data",
    test.use = "wilcox",
    only.pos = TRUE,
    min.pct = MIN_PCT,
    logfc.threshold = LOGFC_THRESHOLD,
    max.cells.per.ident = MAX_CELLS_PER_IDENT_MARKERS,
    random.seed = 20260727,
    verbose = FALSE
  ) %>%
    rownames_to_column("gene")

  lfc_candidates <- c(
    "avg_log2FC",
    "avg_logFC"
  )

  lfc_col <- lfc_candidates[
    lfc_candidates %in% colnames(markers)
  ][1]

  if (is.na(lfc_col)) {
    stop(
      "Cannot detect logFC column in FindMarkers result."
    )
  }

  markers <- markers %>%
    arrange(
      p_val_adj,
      desc(.data[[lfc_col]]),
      desc(pct.1 - pct.2)
    )

  genes <- head(
    markers$gene,
    top_n
  )

  write_csv(
    markers,
    file.path(
      TABLE_DIR,
      "mouse_A4_vs_Other_all_positive_markers_for_signature.csv"
    )
  )

  write_csv(
    tibble(
      rank = seq_along(genes),
      gene = genes
    ),
    file.path(
      TABLE_DIR,
      paste0(
        "mouse_A4_signature_top",
        length(genes),
        ".csv"
      )
    )
  )

  genes
}

p_to_label <- function(p) {
  if (is.na(p)) {
    return("NA")
  }

  if (p < 0.0001) {
    return("****")
  }

  if (p < 0.001) {
    return("***")
  }

  if (p < 0.01) {
    return("**")
  }

  if (p < 0.05) {
    return("*")
  }

  "ns"
}

make_sig_bracket <- function(
  df,
  y_col,
  p_value
) {
  y <- df[[y_col]]
  y <- y[
    is.finite(y)
  ]

  if (length(y) == 0) {
    stop(
      "No finite values found for significance bracket: ",
      y_col
    )
  }

  y_min <- min(
    y,
    na.rm = TRUE
  )

  y_max <- max(
    y,
    na.rm = TRUE
  )

  y_range <- y_max - y_min

  if (
    !is.finite(y_range) ||
    y_range <= 0
  ) {
    y_range <- max(
      abs(y_max),
      0.1
    )
  }

  tibble(
    x1 = 1,
    x2 = 2,
    xm = 1.5,

    y_bracket =
      y_max + 0.08 * y_range,

    y_tick =
      y_max + 0.035 * y_range,

    y_text =
      y_max + 0.145 * y_range,

    y_upper =
      y_max + 0.24 * y_range,

    p_value = p_value,
    label = p_to_label(p_value)
  )
}

# =============================================================================
# 2. 读取对象并检查 metadata
# =============================================================================
if (!file.exists(MOUSE_RDS)) {
  stop(
    "Mouse RDS not found: ",
    MOUSE_RDS
  )
}

cat0("Loading: ", MOUSE_RDS)

obj <- readRDS(MOUSE_RDS)
obj <- join_layers_safe(obj)

if (!"RNA" %in% Assays(obj)) {
  stop("RNA assay not found.")
}

DefaultAssay(obj) <- "RNA"
obj <- add_a4_status(obj)

meta <- obj@meta.data

sample_col <- pick_col(
  meta = meta,
  override = SAMPLE_COL_OVERRIDE,
  candidates = c(
    "sample",
    "original_sample",
    "orig.ident",
    "sample_id",
    "donor"
  ),
  label = "mouse/sample column"
)

genotype_col <- pick_col(
  meta = meta,
  override = GENOTYPE_COL_OVERRIDE,
  candidates = c(
    "genotype",
    "original_genotype",
    "apoe_group",
    "APOE",
    "Genotype"
  ),
  label = "genotype column"
)

obj$mouse_id <- as.character(
  obj@meta.data[[sample_col]]
)

obj$APOE_group <- recode_genotype(
  obj@meta.data[[genotype_col]]
)

cat0("Detected sample column: ", sample_col)
cat0("Detected genotype column: ", genotype_col)

cat0("\nA4 status:")
print(
  table(
    obj$A4_status,
    useNA = "ifany"
  )
)

cat0("\nGenotype:")
print(
  table(
    obj$APOE_group,
    useNA = "ifany"
  )
)

cat0("\nMouse x genotype:")
print(
  table(
    obj$mouse_id,
    obj$APOE_group,
    useNA = "ifany"
  )
)

# 为保证 data layer 可用，重新 LogNormalize
obj <- NormalizeData(
  obj,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

# =============================================================================
# 3. 读取或构建 A4 signature
# =============================================================================
available_genes <- rownames(obj)

if (nzchar(SIGNATURE_FILE)) {
  if (!file.exists(SIGNATURE_FILE)) {
    stop(
      "A4_SIGNATURE_FILE not found: ",
      SIGNATURE_FILE
    )
  }

  signature_genes <- read_signature(
    path = SIGNATURE_FILE,
    available_genes = available_genes
  )

  signature_source <- SIGNATURE_FILE

} else {
  signature_genes <- derive_a4_signature(
    obj = obj,
    top_n = TOP_SIGNATURE_N
  )

  signature_source <-
    "Derived from original mouse A4-OLs vs Other OLs"
}

if (length(signature_genes) < 10) {
  stop(
    "Fewer than 10 usable A4 signature genes: ",
    length(signature_genes)
  )
}

writeLines(
  signature_genes,
  file.path(
    TABLE_DIR,
    "mouse_A4_signature_genes_used.txt"
  )
)

writeLines(
  signature_source,
  file.path(
    TABLE_DIR,
    "mouse_A4_signature_source.txt"
  )
)

cat0(
  "A4 signature genes used: ",
  length(signature_genes)
)

# =============================================================================
# 4. UCell scoring
# =============================================================================
metadata_before_ucell <- colnames(
  obj@meta.data
)

obj <- UCell::AddModuleScore_UCell(
  obj = obj,
  features = list(
    mouse_A4_signature = signature_genes
  ),
  assay = "RNA",
  name = "_UCell"
)

expected_score_col <-
  "mouse_A4_signature_UCell"

if (
  expected_score_col %in%
    colnames(obj@meta.data)
) {
  score_col <- expected_score_col

} else {
  metadata_after_ucell <- colnames(
    obj@meta.data
  )

  new_metadata_cols <- setdiff(
    metadata_after_ucell,
    metadata_before_ucell
  )

  numeric_new_cols <- new_metadata_cols[
    vapply(
      new_metadata_cols,
      function(x) {
        is.numeric(
          obj@meta.data[[x]]
        )
      },
      logical(1)
    )
  ]

  score_candidates <- numeric_new_cols[
    grepl(
      "mouse_A4_signature|A4.*UCell|UCell",
      numeric_new_cols,
      ignore.case = TRUE
    )
  ]

  if (length(score_candidates) == 0) {
    cat0(
      "Metadata columns after UCell:"
    )

    print(
      colnames(obj@meta.data)
    )

    stop(
      "UCell calculation returned without error, ",
      "but the score column could not be identified."
    )
  }

  score_col <- score_candidates[1]
}

cat0(
  "UCell score column: ",
  score_col
)

cat0(
  "UCell score summary:"
)

print(
  summary(
    obj@meta.data[[score_col]]
  )
)

# =============================================================================
# 5. 构建 cell-level 和 mouse-level 数据
# =============================================================================
plot_df <- obj@meta.data %>%
  rownames_to_column("cell") %>%
  transmute(
    cell = cell,
    mouse_id = as.character(mouse_id),

    genotype = factor(
      APOE_group,
      levels = c(
        "APOE3",
        "APOE4"
      )
    ),

    A4_status =
      as.character(A4_status),

    A4_UCell =
      as.numeric(.data[[score_col]])
  ) %>%
  filter(
    A4_status == "A4-OLs",
    genotype %in% c(
      "APOE3",
      "APOE4"
    ),
    is.finite(A4_UCell),
    !is.na(mouse_id),
    mouse_id != ""
  )

if (nrow(plot_df) == 0) {
  stop(
    "No eligible A4-OL cells for APOE3/APOE4 comparison."
  )
}

# -----------------------------------------------------------------------------
# 按每只 mouse 的 A4-OL 细胞数过滤
# 仅保留 n_A4_cells >= 10 的 mouse
# -----------------------------------------------------------------------------
mouse_cell_counts <- plot_df %>%
  count(
    mouse_id,
    genotype,
    name = "n_A4_cells"
  ) %>%
  arrange(
    genotype,
    mouse_id
  )

eligible_mice <- mouse_cell_counts %>%
  filter(
    n_A4_cells >= MIN_A4_CELLS_PER_MOUSE
  )

excluded_mice <- mouse_cell_counts %>%
  filter(
    n_A4_cells < MIN_A4_CELLS_PER_MOUSE
  )

cat0(
  "\nMinimum A4-OL cells per mouse: ",
  MIN_A4_CELLS_PER_MOUSE
)

cat0("\nIncluded mice:")
print(
  as.data.frame(eligible_mice),
  row.names = FALSE
)

cat0("\nExcluded mice:")
print(
  as.data.frame(excluded_mice),
  row.names = FALSE
)

write_csv(
  mouse_cell_counts,
  file.path(
    TABLE_DIR,
    "mouse_A4OL_cell_counts_before_min10_filter.csv"
  )
)

write_csv(
  eligible_mice,
  file.path(
    TABLE_DIR,
    "mouse_A4OL_mice_included_min10.csv"
  )
)

write_csv(
  excluded_mice,
  file.path(
    TABLE_DIR,
    "mouse_A4OL_mice_excluded_below10.csv"
  )
)

# cell-level 和 mouse-level 分析都只使用合格 mouse
plot_df <- plot_df %>%
  semi_join(
    eligible_mice %>%
      select(mouse_id, genotype),
    by = c("mouse_id", "genotype")
  )

if (nrow(plot_df) == 0) {
  stop(
    "No A4-OL cells remain after applying MIN_A4_CELLS_PER_MOUSE = ",
    MIN_A4_CELLS_PER_MOUSE
  )
}

if (n_distinct(plot_df$genotype) < 2) {
  stop(
    "After filtering, both APOE3 and APOE4 groups are required."
  )
}

write_csv(
  plot_df,
  file.path(
    TABLE_DIR,
    "mouse_A4OL_UCell_APOE3_vs_APOE4_per_cell.csv"
  )
)

mouse_df <- plot_df %>%
  group_by(
    mouse_id,
    genotype
  ) %>%
  summarise(
    n_A4_cells = n(),
    mean_A4_UCell = mean(
      A4_UCell,
      na.rm = TRUE
    ),
    median_A4_UCell = median(
      A4_UCell,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

write_csv(
  mouse_df,
  file.path(
    TABLE_DIR,
    "mouse_A4OL_UCell_APOE3_vs_APOE4_per_mouse.csv"
  )
)

cat0("\nPer-mouse A4-OL UCell summary after min10 filtering:")
print(
  mouse_df,
  n = Inf,
  width = Inf
)

# =============================================================================
# 6. 统计检验
# =============================================================================
cell_p <- tryCatch(
  wilcox.test(
    A4_UCell ~ genotype,
    data = plot_df,
    exact = FALSE
  )$p.value,
  error = function(e) {
    cat0(
      "[Warning] Cell-level Wilcoxon failed: ",
      conditionMessage(e)
    )
    NA_real_
  }
)

mouse_p <- tryCatch(
  wilcox.test(
    median_A4_UCell ~ genotype,
    data = mouse_df,
    exact = FALSE
  )$p.value,
  error = function(e) {
    cat0(
      "[Warning] Mouse-level Wilcoxon failed: ",
      conditionMessage(e)
    )
    NA_real_
  }
)

stats <- tibble(
  analysis_level = c(
    "cell_level_descriptive",
    "mouse_level_primary"
  ),

  n_APOE3 = c(
    sum(plot_df$genotype == "APOE3"),
    sum(mouse_df$genotype == "APOE3")
  ),

  n_APOE4 = c(
    sum(plot_df$genotype == "APOE4"),
    sum(mouse_df$genotype == "APOE4")
  ),

  mean_APOE3 = c(
    mean(
      plot_df$A4_UCell[
        plot_df$genotype == "APOE3"
      ],
      na.rm = TRUE
    ),

    mean(
      mouse_df$median_A4_UCell[
        mouse_df$genotype == "APOE3"
      ],
      na.rm = TRUE
    )
  ),

  mean_APOE4 = c(
    mean(
      plot_df$A4_UCell[
        plot_df$genotype == "APOE4"
      ],
      na.rm = TRUE
    ),

    mean(
      mouse_df$median_A4_UCell[
        mouse_df$genotype == "APOE4"
      ],
      na.rm = TRUE
    )
  ),

  delta_APOE4_minus_APOE3 =
    mean_APOE4 - mean_APOE3,

  wilcox_p = c(
    cell_p,
    mouse_p
  ),

  significance = vapply(
    c(cell_p, mouse_p),
    p_to_label,
    character(1)
  )
)

write_csv(
  stats,
  file.path(
    TABLE_DIR,
    "mouse_A4OL_UCell_APOE3_vs_APOE4_stats.csv"
  )
)

cat0("\nStatistical results:")
print(
  stats,
  n = Inf,
  width = Inf
)

# 可选 mixed-effects model
if (
  requireNamespace(
    "lme4",
    quietly = TRUE
  ) &&
  requireNamespace(
    "lmerTest",
    quietly = TRUE
  )
) {
  mixed_fit <- tryCatch(
    lmerTest::lmer(
      A4_UCell ~ genotype + (1 | mouse_id),
      data = plot_df
    ),
    error = function(e) {
      cat0(
        "[Warning] Mixed model failed: ",
        conditionMessage(e)
      )
      NULL
    }
  )

  if (!is.null(mixed_fit)) {
    capture.output(
      summary(mixed_fit),
      file = file.path(
        TABLE_DIR,
        "mouse_A4OL_UCell_mixed_model.txt"
      )
    )
  }
}

# =============================================================================
# 7. 显著性括号数据
# =============================================================================
cell_sig <- make_sig_bracket(
  df = plot_df,
  y_col = "A4_UCell",
  p_value = cell_p
)

mouse_sig <- make_sig_bracket(
  df = mouse_df,
  y_col = "median_A4_UCell",
  p_value = mouse_p
)

write_csv(
  cell_sig,
  file.path(
    TABLE_DIR,
    "mouse_A4OL_UCell_cell_level_bracket.csv"
  )
)

write_csv(
  mouse_sig,
  file.path(
    TABLE_DIR,
    "mouse_A4OL_UCell_mouse_level_bracket.csv"
  )
)

cat0(
  "Cell-level significance label: ",
  cell_sig$label
)

cat0(
  "Mouse-level significance label: ",
  mouse_sig$label
)

# =============================================================================
# 8. 作图
# =============================================================================
base_theme <- theme_classic(
  base_family = PLOT_FONT,
  base_size = 15
) +
  theme(
    text = element_text(
      family = PLOT_FONT,
      color = "black"
    ),

    plot.title = element_text(
      size = TITLE_SIZE,
      face = "bold",
      hjust = 0.5,
      color = "black"
    ),

    plot.subtitle = element_text(
      size = SUBTITLE_SIZE,
      hjust = 0.5,
      color = "black"
    ),

    axis.title = element_text(
      size = AXIS_TITLE_SIZE,
      color = "black"
    ),

    axis.text = element_text(
      size = AXIS_TEXT_SIZE,
      color = "black"
    ),

    legend.position = "none",

    panel.border = element_rect(
      fill = NA,
      color = "black",
      linewidth = 0.7
    ),

    plot.margin = margin(
      t = 12,
      r = 10,
      b = 8,
      l = 8
    )
  )

# -----------------------------------------------------------------------------
# 8.1 Cell-level 图
# -----------------------------------------------------------------------------
p_cell <- ggplot(
  plot_df,
  aes(
    x = genotype,
    y = A4_UCell,
    fill = genotype,
    color = genotype
  )
) +
  geom_violin(
    trim = TRUE,
    alpha = 0.18,
    linewidth = 0.6
  ) +
  geom_jitter(
    width = 0.14,
    size = 0.45,
    alpha = 0.25
  ) +
  geom_boxplot(
    width = 0.20,
    outlier.shape = NA,
    alpha = 0.75,
    linewidth = 0.55
  ) +

  # 显著性横线
  geom_segment(
    data = cell_sig,
    aes(
      x = x1,
      xend = x2,
      y = y_bracket,
      yend = y_bracket
    ),
    inherit.aes = FALSE,
    color = "black",
    linewidth = SIG_LINE_WIDTH
  ) +

  # 左侧竖线
  geom_segment(
    data = cell_sig,
    aes(
      x = x1,
      xend = x1,
      y = y_tick,
      yend = y_bracket
    ),
    inherit.aes = FALSE,
    color = "black",
    linewidth = SIG_LINE_WIDTH
  ) +

  # 右侧竖线
  geom_segment(
    data = cell_sig,
    aes(
      x = x2,
      xend = x2,
      y = y_tick,
      yend = y_bracket
    ),
    inherit.aes = FALSE,
    color = "black",
    linewidth = SIG_LINE_WIDTH
  ) +

  # 星号或 ns
  geom_text(
    data = cell_sig,
    aes(
      x = xm,
      y = y_text,
      label = label
    ),
    inherit.aes = FALSE,
    family = PLOT_FONT,
    color = "black",
    size = SIG_TEXT_SIZE
  ) +

  scale_fill_manual(
    values = COLORS
  ) +

  scale_color_manual(
    values = COLORS
  ) +

  coord_cartesian(
    ylim = c(
      0,
      cell_sig$y_upper
    ),
    clip = "off"
  ) +

  labs(
    title =
      "Mouse A4-OL signature activity by APOE genotype",


    x = NULL,

    y =
      "A4-OL UCell score"
  ) +

  base_theme

# -----------------------------------------------------------------------------
# 8.2 Mouse-level 图
# -----------------------------------------------------------------------------
p_mouse <- ggplot(
  mouse_df,
  aes(
    x = genotype,
    y = median_A4_UCell,
    fill = genotype,
    color = genotype
  )
) +
  geom_boxplot(
    width = 0.40,
    outlier.shape = NA,
    alpha = 0.35,
    linewidth = 0.65
  ) +
  geom_jitter(
    width = 0.08,
    size = 3.0,
    alpha = 0.95
  ) +

  # 显著性横线
  geom_segment(
    data = mouse_sig,
    aes(
      x = x1,
      xend = x2,
      y = y_bracket,
      yend = y_bracket
    ),
    inherit.aes = FALSE,
    color = "black",
    linewidth = SIG_LINE_WIDTH
  ) +

  # 左侧竖线
  geom_segment(
    data = mouse_sig,
    aes(
      x = x1,
      xend = x1,
      y = y_tick,
      yend = y_bracket
    ),
    inherit.aes = FALSE,
    color = "black",
    linewidth = SIG_LINE_WIDTH
  ) +

  # 右侧竖线
  geom_segment(
    data = mouse_sig,
    aes(
      x = x2,
      xend = x2,
      y = y_tick,
      yend = y_bracket
    ),
    inherit.aes = FALSE,
    color = "black",
    linewidth = SIG_LINE_WIDTH
  ) +

  # 星号或 ns
  geom_text(
    data = mouse_sig,
    aes(
      x = xm,
      y = y_text,
      label = label
    ),
    inherit.aes = FALSE,
    family = PLOT_FONT,
    color = "black",
    size = SIG_TEXT_SIZE
  ) +

  scale_fill_manual(
    values = COLORS
  ) +

  scale_color_manual(
    values = COLORS
  ) +

  coord_cartesian(
    ylim = c(
      0,
      mouse_sig$y_upper
    ),
    clip = "off"
  ) +

  labs(
    title =
      "Mouse-level A4-OL signature activity",

    subtitle =
      "Each point represents one independent mouse",

    x = NULL,

    y =
      "Median A4-OL UCell score per mouse"
  ) +

  base_theme

# =============================================================================
# 9. 保存图片
# =============================================================================
ggsave(
  filename = file.path(
    PLOT_DIR,
    "Mouse_A4OL_UCell_APOE3_vs_APOE4_cell_level_with_significance.png"
  ),
  plot = p_cell,
  width = CELL_PLOT_WIDTH,
  height = CELL_PLOT_HEIGHT,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(
    PLOT_DIR,
    "Mouse_A4OL_UCell_APOE3_vs_APOE4_cell_level_with_significance.pdf"
  ),
  plot = p_cell,
  width = CELL_PLOT_WIDTH,
  height = CELL_PLOT_HEIGHT,
  bg = "white"
)

ggsave(
  filename = file.path(
    PLOT_DIR,
    "Mouse_A4OL_UCell_APOE3_vs_APOE4_mouse_level_with_significance.png"
  ),
  plot = p_mouse,
  width = MOUSE_PLOT_WIDTH,
  height = MOUSE_PLOT_HEIGHT,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(
    PLOT_DIR,
    "Mouse_A4OL_UCell_APOE3_vs_APOE4_mouse_level_with_significance.pdf"
  ),
  plot = p_mouse,
  width = MOUSE_PLOT_WIDTH,
  height = MOUSE_PLOT_HEIGHT,
  bg = "white"
)

# 同时保留原文件名，便于覆盖旧图
ggsave(
  filename = file.path(
    PLOT_DIR,
    "Mouse_A4OL_UCell_APOE3_vs_APOE4_cell_level.png"
  ),
  plot = p_cell,
  width = CELL_PLOT_WIDTH,
  height = CELL_PLOT_HEIGHT,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(
    PLOT_DIR,
    "Mouse_A4OL_UCell_APOE3_vs_APOE4_cell_level.pdf"
  ),
  plot = p_cell,
  width = CELL_PLOT_WIDTH,
  height = CELL_PLOT_HEIGHT,
  bg = "white"
)

ggsave(
  filename = file.path(
    PLOT_DIR,
    "Mouse_A4OL_UCell_APOE3_vs_APOE4_mouse_level.png"
  ),
  plot = p_mouse,
  width = MOUSE_PLOT_WIDTH,
  height = MOUSE_PLOT_HEIGHT,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(
    PLOT_DIR,
    "Mouse_A4OL_UCell_APOE3_vs_APOE4_mouse_level.pdf"
  ),
  plot = p_mouse,
  width = MOUSE_PLOT_WIDTH,
  height = MOUSE_PLOT_HEIGHT,
  bg = "white"
)

# =============================================================================
# 10. 保存带 UCell 分数的对象
# =============================================================================
saveRDS(
  obj,
  file.path(
    RDS_DIR,
    "mouse_olig_with_A4_UCell_score.rds"
  ),
  compress = TRUE
)

cat0("\n============================================================")
cat0("DONE")
cat0("============================================================")
cat0("Output: ", OUT_ROOT)

cat0(
  "Cell-level P = ",
  format.pval(
    cell_p,
    digits = 4,
    eps = 1e-4
  ),
  "; label = ",
  p_to_label(cell_p)
)

cat0(
  "Mouse-level P = ",
  format.pval(
    mouse_p,
    digits = 4,
    eps = 1e-4
  ),
  "; label = ",
  p_to_label(mouse_p)
)
