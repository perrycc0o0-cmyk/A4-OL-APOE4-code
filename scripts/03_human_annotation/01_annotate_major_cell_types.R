# ============================================================================
# 手动细胞类型标注（完全修正版）
# 细胞类型名称已改为单数：Excitatory Neuron, Inhibitory Neuron, Astrocyte, Oligodendrocyte
# UMAP 可视化已与 "细胞注释模块-最终版" 统一（无左移）
# 输出到 result-qin/Figure/Mapping/result
#
# 修改：
# 1. celltype 配色改为指定马卡龙配色
# 2. UMAP_cell_type_manual.png 图最上方只标注总细胞数：xxx cells
# 3. 每个 celltype 区域标签仍然只显示 celltype 名称
# ============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(showtext)

source(file.path("R", "load_config.R"))

# ------------------------------ 字体设置 ------------------------------------
font_add("Arial", A4OL_FONT_FILE)
showtext_auto()

theme_set(
  theme_minimal(base_family = "Arial", base_size = 50) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(color = "black", size = 50),
      axis.title = element_text(size = 50),
      plot.title = element_text(size = 50, face = "bold", hjust = 0.5),
      legend.text = element_text(size = 50),
      legend.title = element_text(size = 50)
    )
)

# ------------------------------ 路径设置 ------------------------------------
input_rds <- file.path(
  A4OL_PROJECT_ROOT,
  "20260419",
  "snRNAseq",
  "integrated_harmony_with_celltypes.rds"
)

output_dir <- file.path(A4OL_FIGURE_ROOT, "Mapping", "result")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

# ------------------------------ 标记基因列表（单数名称）----------------------
cell_type_markers <- list(
  "Excitatory Neuron" = c("SATB2", "SLC17A7"),
  "Inhibitory Neuron" = c("GAD1", "GAD2"),
  "Oligodendrocyte"   = c("MOBP", "MAG", "MOG"),
  "OPC"               = c("VCAN", "PDGFRA", "CSPG4", "BCAN"),
  "Astrocyte"         = c("GFAP", "AQP4", "GJA1", "ALDH1L1"),
  "Microglia"         = c("CSF1R", "CD74", "C3"),
  "Vascular"          = c("CLDN5", "FLT1", "PDGFRB")
)

# -------------------------- 手动 cluster 映射表（单数名称）------------------
cluster_to_celltype <- list(
  "Excitatory Neuron" = c(8, 33, 19, 14, 34, 20, 22, 28, 9, 6, 18, 26, 13),
  "Inhibitory Neuron" = c(10, 21, 15, 23, 36, 30),
  "Oligodendrocyte"   = c(25, 7, 0, 42, 32),
  "OPC"               = c(44, 29, 11),
  "Astrocyte"         = c(12, 46, 1, 4, 38, 5, 27, 24),
  "Microglia"         = c(2, 3, 45, 35, 41, 31),
  "Vascular"          = c(16, 39, 37, 17, 40, 43)
)

# ------------------------------ 颜色方案：指定 celltype 配色 ----------------
enhanced_cell_type_colors <- c(
  "Excitatory Neuron"    = "#F7A6AC",  # 粉红
  "Inhibitory Neuron"    = "#EEC78A",  # 杏黄
  "Oligodendrocyte"      = "#B8E5FA",  # 浅蓝
  "OPC"                  = "#B3DDCB",  # 青绿
  "Astrocyte"            = "#CBE4B1",  # 浅绿
  "Microglia"            = "#EEE9A2",  # 浅黄
  "Vascular"             = "#F7B2C7",  # 玫粉
  "Unknown"              = "#D9D9D9"   # 灰色
)

# ------------------------------ 1. 加载数据 --------------------------------
cat("加载 RDS 文件...\n")

scRNA <- readRDS(input_rds)
DefaultAssay(scRNA) <- "RNA"

cat("细胞总数:", ncol(scRNA), "\n")

if (!"seurat_clusters" %in% colnames(scRNA@meta.data)) {
  stop("RDS 中没有 seurat_clusters 列！")
}

actual_clusters <- sort(unique(scRNA$seurat_clusters))

cat("实际存在的 seurat_clusters 编号:", paste(actual_clusters, collapse = ", "), "\n")

# ------------------------------ 2. 应用手动标注 -----------------------------
mapping_df <- stack(cluster_to_celltype)
colnames(mapping_df) <- c("cluster", "cell_type")
mapping_df$cluster <- as.character(mapping_df$cluster)

scRNA$cell_type_manual <- NA

for (i in seq_len(nrow(mapping_df))) {
  clus <- mapping_df$cluster[i]
  ct   <- as.character(mapping_df$cell_type[i])
  scRNA$cell_type_manual[scRNA$seurat_clusters == clus] <- ct
}

# 处理未映射的 cluster
unmapped <- setdiff(actual_clusters, mapping_df$cluster)

if (length(unmapped) > 0) {
  cat("未映射的 cluster:", paste(unmapped, collapse = ", "), " -> 标记为 Unknown\n")
  for (clus in unmapped) {
    scRNA$cell_type_manual[scRNA$seurat_clusters == clus] <- "Unknown"
  }
}

# 设置因子水平
cell_type_order <- names(cell_type_markers)

if ("Unknown" %in% unique(scRNA$cell_type_manual)) {
  cell_type_order <- c(cell_type_order, "Unknown")
}

scRNA$cell_type_manual <- factor(scRNA$cell_type_manual, levels = cell_type_order)

cat("\n手动标注统计:\n")
print(table(scRNA$cell_type_manual))

# ------------------------------ 3. 可视化 -----------------------------------
main_reduction <- if ("umap.harmony" %in% names(scRNA@reductions)) {
  "umap.harmony"
} else {
  "umap"
}

# 构建命名颜色向量
plot_colors <- enhanced_cell_type_colors[levels(scRNA$cell_type_manual)]
plot_colors <- plot_colors[!is.na(plot_colors)]

# ---------- UMAP 箭头函数（字体 50，UMAP-2 上移 *8）----------
add_umap_arrows <- function(obj, p, reduction_name) {
  emb <- Embeddings(obj, reduction_name)

  x_range <- range(emb[, 1])
  y_range <- range(emb[, 2])

  arrow_start_x <- x_range[1] - diff(x_range) * 0.05
  arrow_start_y <- y_range[1] - diff(y_range) * 0.05

  len_x <- diff(x_range) * 0.15
  len_y <- diff(y_range) * 0.15

  offset_x <- diff(x_range) * 0.03
  offset_y <- diff(y_range) * 0.03

  p +
    theme(
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    ) +
    annotate(
      "segment",
      x = arrow_start_x,
      xend = arrow_start_x + len_x,
      y = arrow_start_y,
      yend = arrow_start_y,
      arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
      linewidth = 1.2,
      color = "black"
    ) +
    annotate(
      "segment",
      x = arrow_start_x,
      xend = arrow_start_x,
      y = arrow_start_y,
      yend = arrow_start_y + len_y,
      arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
      linewidth = 1.2,
      color = "black"
    ) +
    annotate(
      "text",
      x = arrow_start_x + offset_x,
      y = arrow_start_y - offset_y,
      label = "UMAP-1",
      hjust = 0,
      vjust = 1,
      size = 50,
      family = "Arial"
    ) +
    annotate(
      "text",
      x = arrow_start_x - offset_x,
      y = arrow_start_y + offset_y * 8,
      label = "UMAP-2",
      hjust = 1,
      vjust = 0,
      angle = 90,
      size = 50,
      family = "Arial"
    ) +
    coord_cartesian(clip = "off")
}

# ============ 自定义 UMAP 细胞类型图（无左移，顶部标注总细胞数）============
umap_emb <- as.data.frame(Embeddings(scRNA, main_reduction))
colnames(umap_emb) <- c("UMAP1", "UMAP2")
umap_emb$cell_type <- scRNA$cell_type_manual

centers <- umap_emb %>%
  group_by(cell_type) %>%
  summarise(
    x = median(UMAP1),
    y = median(UMAP2),
    .groups = "drop"
  )

# 整张图最上方标注总细胞数
total_cells <- ncol(scRNA)

x_mid <- mean(range(umap_emb$UMAP1))
y_top <- max(umap_emb$UMAP2) + diff(range(umap_emb$UMAP2)) * 0.08

p1 <- ggplot(umap_emb, aes(x = UMAP1, y = UMAP2, color = cell_type)) +
  geom_point(size = 0.3, alpha = 0.6) +
  scale_color_manual(values = plot_colors, drop = FALSE) +
  geom_text(
    data = centers,
    aes(x = x, y = y, label = cell_type),
    color = "black",
    size = 120 / ggplot2::.pt,
    family = "Arial",
    fontface = "bold",
    lineheight = 0.9,
    inherit.aes = FALSE
  ) +
  annotate(
    "text",
    x = x_mid,
    y = y_top,
    label = paste0(total_cells, " cells"),
    color = "black",
    size = 120 / ggplot2::.pt,
    family = "Arial",
    fontface = "bold"
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    plot.margin = margin(t = 30, r = 10, b = 10, l = 10)
  )

p1 <- add_umap_arrows(scRNA, p1, main_reduction)

ggsave(
  file.path(output_dir, "plots", "UMAP_cell_type_manual.png"),
  p1,
  width = 12,
  height = 12,
  dpi = 300,
  bg = "white"
)

# UMAP：原始聚类编号（保留以供参考）
p2 <- DimPlot(
  scRNA,
  reduction = main_reduction,
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.3
) +
  ggtitle("Clusters (res = 0.8)")

ggsave(
  file.path(output_dir, "plots", "UMAP_clusters.png"),
  p2,
  width = 14,
  height = 12,
  dpi = 300,
  bg = "white"
)

# DotPlot：红色、无标题、无坐标轴标题
all_markers <- unique(unlist(cell_type_markers))
available_markers <- all_markers[all_markers %in% rownames(scRNA)]
missing_markers <- setdiff(all_markers, available_markers)

if (length(missing_markers) > 0) {
  cat("缺失的标记基因:", paste(missing_markers, collapse = ", "), "\n")
}

if (length(available_markers) > 0) {
  p3 <- DotPlot(
    scRNA,
    features = available_markers,
    group.by = "cell_type_manual",
    cols = c("lightgrey", "#DC2626"),
    dot.scale = 8,
    scale = TRUE
  ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 70, face = "bold"),
      axis.text.y = element_text(size = 70, face = "bold"),
      legend.text = element_text(size = 70),
      legend.title = element_text(size = 70, face = "bold"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_blank()
    ) +
    labs(x = NULL, y = NULL)

  ggsave(
    file.path(output_dir, "plots", "dotplot_all_markers.png"),
    p3,
    width = 16,
    height = 7,
    dpi = 300,
    bg = "white"
  )
}

# 细胞组成比例图（若存在 apoe_group）
if ("apoe_group" %in% colnames(scRNA@meta.data)) {
  cell_comp_sample <- scRNA@meta.data %>%
    group_by(orig.ident, apoe_group, cell_type_manual) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(orig.ident) %>%
    mutate(percentage = count / sum(count) * 100)

  p4 <- ggplot(
    cell_comp_sample,
    aes(x = orig.ident, y = percentage, fill = cell_type_manual)
  ) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid(. ~ apoe_group, scales = "free_x", space = "free") +
    scale_fill_manual(values = plot_colors) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "Cell Type Composition by Sample (Manual Annotation)",
      x = "Sample",
      y = "Percentage (%)"
    )

  ggsave(
    file.path(output_dir, "plots", "cell_type_composition_by_sample_manual.png"),
    p4,
    width = 16,
    height = 8,
    dpi = 300,
    bg = "white"
  )
} else {
  cat("元数据中无 'apoe_group' 列，跳过样本组成比例图。\n")
}

# ------------------------------ 4. 保存结果 --------------------------------
saveRDS(
  scRNA,
  file = file.path(output_dir, "integrated_harmony_manual_celltypes.rds"),
  compress = TRUE
)

write.csv(
  scRNA@meta.data,
  file = file.path(output_dir, "tables", "cell_metadata_manual_annotation.csv")
)

write.csv(
  mapping_df,
  file = file.path(output_dir, "tables", "cluster_to_celltype_mapping.csv"),
  row.names = FALSE
)

cat("\n全部完成！结果保存至:", output_dir, "\n")

showtext_auto(FALSE)
