#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

options(stringsAsFactors = FALSE)
options(repos = c(CRAN = "https://cloud.r-project.org"))

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(readr)
  library(stringr)
})

source(file.path("R", "load_config.R"))

# ============================================================
# 0. Paths
# ============================================================
BASE_DIR <- file.path(A4OL_FIGURE_ROOT, "UCell_A4OL")

UCell_RDS <- file.path(
  BASE_DIR,
  "rds",
  "human_olig_res04_with_A4_UCell.rds"
)

OUT_DIR <- file.path(BASE_DIR, "plots_paper_style_significance_final")
TABLE_DIR <- file.path(BASE_DIR, "tables_paper_style_significance_final")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. Parameters
# ============================================================
set.seed(123)

# 散点只用于展示抽样；统计仍然使用全部细胞
MAX_POINTS_PER_GROUP_FOR_DOTS <- 3000

# 文章风格颜色
GENOTYPE_COLORS <- c(
  "APOE3" = "#4DBBD5",
  "APOE4" = "#E64B35"
)

A4_GROUP_COLORS <- c(
  "Other OLs" = "#71A682",
  "A4-OLs"   = "#D19246"
)

# 图层透明度
POINT_ALPHA  <- 0.35
BOX_ALPHA    <- 0.90
VIOLIN_ALPHA <- 0.50

# ============================================================
# 2. Helper functions
# ============================================================
pick_first_existing <- function(candidates, meta) {
  hit <- candidates[candidates %in% colnames(meta)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

detect_score_col <- function(meta) {
  candidates <- c(
    "A4_OL_UCell",
    "A4_OLs_UCell",
    "A4OL_UCell",
    "A4_UCell",
    "A4_OL_signature_UCell"
  )

  hit <- candidates[candidates %in% colnames(meta)]
  if (length(hit) > 0) return(hit[1])

  hit2 <- grep("A4.*UCell|UCell.*A4", colnames(meta), value = TRUE)
  if (length(hit2) > 0) return(hit2[1])

  return(NA_character_)
}

recode_apoe <- function(x) {
  x <- as.character(x)
  x_upper <- toupper(x)

  out <- ifelse(
    grepl("4", x_upper) | grepl("APOE4", x_upper) | grepl("E4", x_upper),
    "APOE4",
    ifelse(
      grepl("3", x_upper) | grepl("APOE3", x_upper) | grepl("E3", x_upper),
      "APOE3",
      NA_character_
    )
  )

  out
}

p_to_star <- function(p) {
  if (is.na(p)) return("ns")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return("ns")
}

safe_wilcox_p <- function(df, value_col, group_col) {
  df <- as.data.frame(df)
  df <- df[!is.na(df[[value_col]]) & !is.na(df[[group_col]]), , drop = FALSE]

  if (nrow(df) < 2) return(NA_real_)

  group_values <- unique(as.character(df[[group_col]]))
  group_values <- group_values[!is.na(group_values)]

  if (length(group_values) != 2) return(NA_real_)

  tb <- table(df[[group_col]])
  if (any(tb < 1)) return(NA_real_)

  tryCatch({
    wilcox.test(df[[value_col]] ~ df[[group_col]])$p.value
  }, error = function(e) {
    NA_real_
  })
}

sample_points_for_plot <- function(df, group_col, max_n = 3000) {
  df <- as.data.frame(df)

  split_list <- split(df, df[[group_col]])

  sampled_list <- lapply(split_list, function(x) {
    if (nrow(x) > max_n) {
      x[sample(seq_len(nrow(x)), max_n), , drop = FALSE]
    } else {
      x
    }
  })

  bind_rows(sampled_list)
}

make_bracket_df <- function(df, x1, x2, y_col, p_value, y_expand = 0.11) {
  y_max <- max(df[[y_col]], na.rm = TRUE)
  y_min <- min(df[[y_col]], na.rm = TRUE)
  y_range <- y_max - y_min

  if (!is.finite(y_range) || y_range <= 0) {
    y_range <- 1
  }

  data.frame(
    x1 = x1,
    x2 = x2,
    y = y_max + y_range * y_expand,
    y0 = y_max + y_range * (y_expand - 0.035),
    label = p_to_star(p_value),
    stringsAsFactors = FALSE
  )
}

theme_paper <- function(base_size = 15) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(
        angle = 0,
        hjust = 0.5,
        vjust = 0.5,
        face = "bold"
      ),
      axis.title = element_text(color = "black", face = "bold"),
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        color = "black"
      ),
      legend.position = "none",
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.6
      ),
      axis.line = element_line(
        color = "black",
        linewidth = 0.5
      )
    )
}

add_significance_bracket <- function(p, bracket_df, y_range) {
  if (!is.finite(y_range) || y_range <= 0) {
    y_range <- 1
  }

  p +
    geom_segment(
      data = bracket_df,
      aes(x = x1, xend = x2, y = y, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.60,
      color = "black"
    ) +
    geom_segment(
      data = bracket_df,
      aes(x = x1, xend = x1, y = y0, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.60,
      color = "black"
    ) +
    geom_segment(
      data = bracket_df,
      aes(x = x2, xend = x2, y = y0, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.60,
      color = "black"
    ) +
    geom_text(
      data = bracket_df,
      aes(x = 1.5, y = y + y_range * 0.018, label = label),
      inherit.aes = FALSE,
      size = 6,
      fontface = "bold",
      color = "black"
    )
}

# ============================================================
# 3. Load UCell object
# ============================================================
cat("Loading UCell RDS:\n", UCell_RDS, "\n")

if (!file.exists(UCell_RDS)) {
  stop("UCell RDS not found: ", UCell_RDS)
}

obj <- readRDS(UCell_RDS)
meta <- obj@meta.data

score_col <- detect_score_col(meta)

cluster_col <- pick_first_existing(
  c(
    "olig_anno_cluster",
    "annotated_cluster",
    "cluster_annotation",
    "manual_cluster",
    "celltype",
    "cell_type"
  ),
  meta
)

genotype_col <- pick_first_existing(
  c(
    "genotype_use",
    "apoe_group",
    "APOE",
    "genotype",
    "Genotype"
  ),
  meta
)

if (is.na(score_col)) {
  stop("Cannot detect A4 UCell score column.")
}

if (is.na(cluster_col)) {
  stop("Cannot detect cluster annotation column.")
}

if (is.na(genotype_col)) {
  stop("Cannot detect APOE genotype column.")
}

cat("Detected score column   :", score_col, "\n")
cat("Detected cluster column :", cluster_col, "\n")
cat("Detected genotype column:", genotype_col, "\n")

# ============================================================
# 4. Build plotting dataframe
# ============================================================
plot_df <- meta %>%
  rownames_to_column("cell") %>%
  transmute(
    cell = cell,
    cluster = as.character(.data[[cluster_col]]),
    genotype_raw = as.character(.data[[genotype_col]]),
    genotype = recode_apoe(.data[[genotype_col]]),
    A4_UCell = as.numeric(.data[[score_col]])
  ) %>%
  filter(!is.na(cluster), !is.na(A4_UCell))

write_csv(
  plot_df,
  file.path(TABLE_DIR, "all_cell_plot_table_A4_UCell.csv")
)

cat("\nCluster table:\n")
print(table(plot_df$cluster))

cat("\nGenotype table:\n")
print(table(plot_df$genotype, useNA = "ifany"))

# ============================================================
# 5. Figure 1:
#    human_A4_like_OL_A4_UCell_APOE3_vs_APOE4
# ============================================================
cat("\nPlotting Figure 1: A4-OLs APOE3 vs APOE4...\n")

df_a4_apoe <- plot_df %>%
  filter(cluster == "A4-OLs") %>%
  filter(genotype %in% c("APOE3", "APOE4")) %>%
  mutate(
    genotype = factor(genotype, levels = c("APOE3", "APOE4"))
  )

df_a4_apoe_dots <- sample_points_for_plot(
  df = df_a4_apoe,
  group_col = "genotype",
  max_n = MAX_POINTS_PER_GROUP_FOR_DOTS
)

p_a4_apoe <- safe_wilcox_p(
  df = df_a4_apoe,
  value_col = "A4_UCell",
  group_col = "genotype"
)

stat_a4_apoe <- data.frame(
  comparison = "A4-OLs_APOE3_vs_APOE4",
  p_value = p_a4_apoe,
  significance = p_to_star(p_a4_apoe),
  n_APOE3 = sum(df_a4_apoe$genotype == "APOE3"),
  n_APOE4 = sum(df_a4_apoe$genotype == "APOE4"),
  stringsAsFactors = FALSE
)

write_csv(
  df_a4_apoe,
  file.path(TABLE_DIR, "table_A4OLs_APOE3_vs_APOE4_cell_level.csv")
)

write_csv(
  stat_a4_apoe,
  file.path(TABLE_DIR, "stats_A4OLs_APOE3_vs_APOE4_cell_level.csv")
)

bracket_a4_apoe <- make_bracket_df(
  df = df_a4_apoe,
  x1 = "APOE3",
  x2 = "APOE4",
  y_col = "A4_UCell",
  p_value = p_a4_apoe,
  y_expand = 0.11
)

yrange1 <- diff(range(df_a4_apoe$A4_UCell, na.rm = TRUE))
if (!is.finite(yrange1) || yrange1 <= 0) {
  yrange1 <- 1
}

fig1 <- ggplot(
  df_a4_apoe,
  aes(x = genotype, y = A4_UCell, fill = genotype)
) +
  # 第一层：细胞散点
  geom_jitter(
    data = df_a4_apoe_dots,
    aes(x = genotype, y = A4_UCell, color = genotype),
    width = 0.16,
    size = 0.55,
    alpha = POINT_ALPHA,
    inherit.aes = FALSE
  ) +
  # 第二层：箱线图，放在散点之后，避免被散点盖住
  geom_boxplot(
    width = 0.22,
    outlier.shape = NA,
    alpha = BOX_ALPHA,
    color = "black",
    linewidth = 0.65
  ) +
  # 第三层：小提琴图，最后画，低透明度作为分布背景
  geom_violin(
    width = 0.88,
    trim = FALSE,
    alpha = VIOLIN_ALPHA,
    color = "black",
    linewidth = 0.35,
    scale = "width"
  ) +
  scale_fill_manual(values = GENOTYPE_COLORS) +
  scale_color_manual(values = GENOTYPE_COLORS) +
  coord_cartesian(
    ylim = c(
      min(df_a4_apoe$A4_UCell, na.rm = TRUE),
      max(bracket_a4_apoe$y, na.rm = TRUE) + yrange1 * 0.08
    )
  ) +
  labs(
    title = "A4-OLs",
    x = NULL,
    y = "A4-OL UCell score"
  ) +
  theme_paper(base_size = 15)

fig1 <- add_significance_bracket(
  p = fig1,
  bracket_df = bracket_a4_apoe,
  y_range = yrange1
)

ggsave(
  file.path(
    OUT_DIR,
    "human_A4_like_OL_A4_UCell_APOE3_vs_APOE4_paper_style_sig_final.png"
  ),
  fig1,
  width = 5.2,
  height = 5.8,
  dpi = 300
)

ggsave(
  file.path(
    OUT_DIR,
    "human_A4_like_OL_A4_UCell_APOE3_vs_APOE4_paper_style_sig_final.pdf"
  ),
  fig1,
  width = 5.2,
  height = 5.8
)

# ============================================================
# 6. Figure 2:
#    human_A4_UCell_A4OLs_vs_OtherOLs
# ============================================================
cat("\nPlotting Figure 2: A4-OLs vs Other OLs...\n")

df_a4_other <- plot_df %>%
  mutate(
    A4_group = ifelse(cluster == "A4-OLs", "A4-OLs", "Other OLs"),
    A4_group = factor(A4_group, levels = c("Other OLs", "A4-OLs"))
  ) %>%
  filter(!is.na(A4_group), !is.na(A4_UCell))

df_a4_other_dots <- sample_points_for_plot(
  df = df_a4_other,
  group_col = "A4_group",
  max_n = MAX_POINTS_PER_GROUP_FOR_DOTS
)

p_a4_other <- safe_wilcox_p(
  df = df_a4_other,
  value_col = "A4_UCell",
  group_col = "A4_group"
)

stat_a4_other <- data.frame(
  comparison = "A4-OLs_vs_OtherOLs",
  p_value = p_a4_other,
  significance = p_to_star(p_a4_other),
  n_Other_OLs = sum(df_a4_other$A4_group == "Other OLs"),
  n_A4_OLs = sum(df_a4_other$A4_group == "A4-OLs"),
  stringsAsFactors = FALSE
)

write_csv(
  df_a4_other,
  file.path(TABLE_DIR, "table_A4OLs_vs_OtherOLs_cell_level.csv")
)

write_csv(
  stat_a4_other,
  file.path(TABLE_DIR, "stats_A4OLs_vs_OtherOLs_cell_level.csv")
)

bracket_a4_other <- make_bracket_df(
  df = df_a4_other,
  x1 = "Other OLs",
  x2 = "A4-OLs",
  y_col = "A4_UCell",
  p_value = p_a4_other,
  y_expand = 0.11
)

yrange2 <- diff(range(df_a4_other$A4_UCell, na.rm = TRUE))
if (!is.finite(yrange2) || yrange2 <= 0) {
  yrange2 <- 1
}

fig2 <- ggplot(
  df_a4_other,
  aes(x = A4_group, y = A4_UCell, fill = A4_group)
) +
  # 第一层：细胞散点
  geom_jitter(
    data = df_a4_other_dots,
    aes(x = A4_group, y = A4_UCell, color = A4_group),
    width = 0.16,
    size = 0.50,
    alpha = POINT_ALPHA,
    inherit.aes = FALSE
  ) +
  # 第二层：箱线图，放在散点之后，避免被散点盖住
  geom_boxplot(
    width = 0.22,
    outlier.shape = NA,
    alpha = BOX_ALPHA,
    color = "black",
    linewidth = 0.65
  ) +
  # 第三层：小提琴图，最后画，低透明度作为分布背景
  geom_violin(
    width = 0.88,
    trim = FALSE,
    alpha = VIOLIN_ALPHA,
    color = "black",
    linewidth = 0.35,
    scale = "width"
  ) +
  scale_fill_manual(values = A4_GROUP_COLORS) +
  scale_color_manual(values = A4_GROUP_COLORS) +
  coord_cartesian(
    ylim = c(
      min(df_a4_other$A4_UCell, na.rm = TRUE),
      max(bracket_a4_other$y, na.rm = TRUE) + yrange2 * 0.08
    )
  ) +
  labs(
    title = "A4-OLs vs Other OLs",
    x = NULL,
    y = "A4-OL UCell score"
  ) +
  theme_paper(base_size = 15)

fig2 <- add_significance_bracket(
  p = fig2,
  bracket_df = bracket_a4_other,
  y_range = yrange2
)

ggsave(
  file.path(
    OUT_DIR,
    "human_A4_UCell_A4OLs_vs_OtherOLs_paper_style_sig_final.png"
  ),
  fig2,
  width = 5.2,
  height = 5.8,
  dpi = 300
)

ggsave(
  file.path(
    OUT_DIR,
    "human_A4_UCell_A4OLs_vs_OtherOLs_paper_style_sig_final.pdf"
  ),
  fig2,
  width = 5.2,
  height = 5.8
)

# ============================================================
# 7. Done
# ============================================================
cat("\nDone.\n")
cat("Plots saved to:\n", OUT_DIR, "\n")
cat("Tables saved to:\n", TABLE_DIR, "\n")

cat("\nGenerated plot files:\n")
cat(file.path(
  OUT_DIR,
  "human_A4_like_OL_A4_UCell_APOE3_vs_APOE4_paper_style_sig_final.png"
), "\n")

cat(file.path(
  OUT_DIR,
  "human_A4_UCell_A4OLs_vs_OtherOLs_paper_style_sig_final.png"
), "\n")
