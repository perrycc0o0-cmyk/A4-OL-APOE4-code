# ============================================================================
# Oligodendrocytes 亚群分析 + 合并注释 cluster + A4-OLs vs Other DEG/富集
# 输出目录：Figure/olig-c
#
# 主要修改：
# 1. cluster 1 和 4 合并为 A4-OLs
# 2. cluster 0、2、3 分别标注为 Oligo1、Oligo2、Oligo3
# 3. UMAP / sample 分布图 / genotype 分布图全部使用标注后的 cluster
# 4. subcluster 颜色使用“柔绿森林”配色
# 5. 火山图改为 A4-OLs vs Other
# 6. DEG 前不按基因表达比例和 logFC 过滤：
#    min.pct = 0, logfc.threshold = 0, return.thresh = Inf
# 7. 火山图标注上调前 3、下调前 3，并额外标注 Egr1
# ============================================================================

options(timeout = 180)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(RColorBrewer)
  library(showtext)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(ggrepel)
  library(purrr)
})

source(file.path("R", "load_config.R"))

# ---------- 字体设置 ----------
font_path <- A4OL_FONT_FILE

if (file.exists(font_path)) {
  font_add("Arial", font_path)
  showtext_auto()
} else {
  warning("Arial.ttf not found. Use default font instead.")
}

theme_set(
  theme_minimal(base_family = "Arial", base_size = 50) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 50),
      legend.text = element_text(family = "Arial", size = 50),
      legend.title = element_text(family = "Arial", size = 50)
    )
)

# ---------- 参数 ----------
input_rds <- file.path(
  A4OL_RESULTS_ROOT,
  "Annotation-final",
  "seurat_with_cell_type_annotation.rds"
)

output_dir <- file.path(A4OL_FIGURE_ROOT, "olig-c")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

high_resolution <- 0.15
n_pcs <- 30
set.seed(42)

# ---------------------------- Subcluster 配色：柔绿森林 ----------------------------
subcluster_palette_olig <- c(
  "#B8DBB3",
  "#86BC79",
  "#71A682",
  "#81989B",
  "#D19246",
  "#B5AF8B",
  "#7EA4B6",
  "#4A4F7E"
)

# 标注后的 olig cluster 颜色
# A4-OLs 用柔绿森林里的橙色，便于和 Oligo1/2/3 区分
oligo_anno_colors <- c(
  "Oligo1" = "#B8DBB3",
  "Oligo2" = "#86BC79",
  "Oligo3" = "#71A682",
  "A4-OLs" = "#D19246"
)

annotated_levels <- c("Oligo1", "Oligo2", "Oligo3", "A4-OLs")

# ---------- 加载数据 ----------
cat("Loading data...\n")
scRNA <- readRDS(input_rds)
cat(sprintf("Total cells: %d, Genes: %d\n", ncol(scRNA), nrow(scRNA)))

if (!"cell_type" %in% colnames(scRNA@meta.data)) {
  stop("No 'cell_type' column found in scRNA@meta.data.")
}

# 兼容 Oligodendrocytes / Oligodendrocyte 两种命名
target_cell_type_candidates <- c("Oligodendrocytes", "Oligodendrocyte")
target_cell_type <- target_cell_type_candidates[
  target_cell_type_candidates %in% unique(as.character(scRNA$cell_type))
]

if (length(target_cell_type) == 0) {
  cat("Available cell_type values:\n")
  print(unique(as.character(scRNA$cell_type)))
  stop("Cannot find Oligodendrocytes / Oligodendrocyte in cell_type.")
}

target_cell_type <- target_cell_type[1]
cat("Target cell type:", target_cell_type, "\n")

oligo_cells <- colnames(scRNA)[as.character(scRNA$cell_type) == target_cell_type]

if (length(oligo_cells) == 0) {
  stop("No oligodendrocyte cells found.")
}

oligo_data <- subset(scRNA, cells = oligo_cells)
cat(sprintf("Extracted %d %s cells\n", ncol(oligo_data), target_cell_type))

if ("sample" %in% colnames(oligo_data@meta.data)) {
  oligo_data$original_sample <- oligo_data$sample
}

if ("genotype" %in% colnames(oligo_data@meta.data)) {
  oligo_data$original_genotype <- oligo_data$genotype
}

if ("RNA" %in% Assays(oligo_data)) {
  DefaultAssay(oligo_data) <- "RNA"
}

# ---------- 重新处理 ----------
cat("\nPreprocessing oligodendrocytes...\n")

oligo_data <- NormalizeData(
  oligo_data,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

oligo_data <- FindVariableFeatures(
  oligo_data,
  selection.method = "vst",
  nfeatures = 3000,
  verbose = FALSE
)

oligo_data <- ScaleData(
  oligo_data,
  features = rownames(oligo_data),
  verbose = FALSE
)

oligo_data <- RunPCA(
  oligo_data,
  features = VariableFeatures(oligo_data),
  npcs = 50,
  verbose = FALSE
)

elbow_plot <- ElbowPlot(oligo_data, ndims = 50) +
  ggtitle("PCA Elbow Plot (Oligodendrocytes)")

ggsave(
  file.path(output_dir, "plots/elbow_plot.png"),
  elbow_plot,
  width = 6,
  height = 5,
  dpi = 300,
  bg = "white"
)

# ---------- Harmony ----------
if ("sample" %in% colnames(oligo_data@meta.data)) {
  if (requireNamespace("harmony", quietly = TRUE)) {
    library(harmony)
    oligo_data <- RunHarmony(
      oligo_data,
      group.by.vars = "sample",
      dims.use = 1:n_pcs,
      verbose = TRUE
    )
    cat("Harmony finished.\n")
  } else {
    warning("Package 'harmony' not found. Use PCA directly.")
  }
}

reduction_use <- if ("harmony" %in% names(oligo_data@reductions)) {
  "harmony"
} else {
  "pca"
}

reduction_name <- if ("harmony" %in% names(oligo_data@reductions)) {
  "umap.harmony"
} else {
  "umap.pca"
}

# ---------- UMAP / 聚类 ----------
cat("\nRunning UMAP and clustering...\n")

oligo_data <- RunUMAP(
  oligo_data,
  reduction = reduction_use,
  dims = 1:n_pcs,
  reduction.name = reduction_name,
  reduction.key = "UMAP_",
  verbose = FALSE
)

oligo_data <- FindNeighbors(
  oligo_data,
  reduction = reduction_use,
  dims = 1:n_pcs,
  verbose = FALSE
)

oligo_data <- FindClusters(
  oligo_data,
  resolution = high_resolution,
  verbose = FALSE,
  random.seed = 42
)

cluster_column <- paste0("RNA_snn_res.", high_resolution)

if (!cluster_column %in% colnames(oligo_data@meta.data)) {
  stop("Cannot find clustering column: ", cluster_column)
}

oligo_data$oligo_clusters <- as.character(oligo_data[[cluster_column]][, 1])

sort_cluster_ids <- function(x) {
  x <- unique(as.character(x))
  if (all(grepl("^[0-9]+$", x))) {
    as.character(sort(as.numeric(x)))
  } else {
    sort(x)
  }
}

cluster_levels <- sort_cluster_ids(oligo_data$oligo_clusters)
n_clusters <- length(cluster_levels)

cat(sprintf("Clusters identified: %d\n", n_clusters))
cat("Original clusters:", paste(cluster_levels, collapse = ", "), "\n")

# ---------- cluster 合并注释 ----------
# cluster 0 -> Oligo1
# cluster 1 -> A4-OLs
# cluster 2 -> Oligo2
# cluster 3 -> Oligo3
# cluster 4 -> A4-OLs
cluster_annotation_map <- c(
  "0" = "Oligo1",
  "1" = "A4-OLs",
  "2" = "Oligo2",
  "3" = "Oligo3",
  "4" = "A4-OLs"
)

unmapped_clusters <- setdiff(cluster_levels, names(cluster_annotation_map))

if (length(unmapped_clusters) > 0) {
  warning(
    "These clusters are not in cluster_annotation_map: ",
    paste(unmapped_clusters, collapse = ", "),
    ". They will be named as ClusterX."
  )
  for (cl in unmapped_clusters) {
    cluster_annotation_map[cl] <- paste0("Cluster", cl)
  }
}

oligo_data$oligo_annotated_cluster <- unname(
  cluster_annotation_map[as.character(oligo_data$oligo_clusters)]
)

extra_anno <- setdiff(unique(oligo_data$oligo_annotated_cluster), annotated_levels)
if (length(extra_anno) > 0) {
  extra_colors <- setNames(
    colorRampPalette(subcluster_palette_olig)(length(extra_anno)),
    extra_anno
  )
  oligo_anno_colors <- c(oligo_anno_colors, extra_colors)
  annotated_levels <- c(annotated_levels, extra_anno)
}

oligo_data$oligo_annotated_cluster <- factor(
  oligo_data$oligo_annotated_cluster,
  levels = annotated_levels
)

# 保存 cluster 注释映射表
cluster_annotation_table <- data.frame(
  original_cluster = names(cluster_annotation_map),
  annotated_cluster = unname(cluster_annotation_map),
  stringsAsFactors = FALSE
) %>%
  arrange(
    suppressWarnings(as.numeric(original_cluster))
  )

write.csv(
  cluster_annotation_table,
  file.path(output_dir, "tables/oligo_cluster_annotation_mapping.csv"),
  row.names = FALSE
)

# ---------- UMAP 箭头函数：沿用原图格式 ----------
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
      size = 30,
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
      size = 30,
      family = "Arial"
    ) +
    coord_cartesian(clip = "off")
}

# ============================================================================
# 1. UMAP：标注后的 olig cluster
# 文件名沿用 UMAP_oligo_clusters.png，但保存到 olig-c
# ============================================================================

p_anno_umap <- DimPlot(
  oligo_data,
  reduction = reduction_name,
  group.by = "oligo_annotated_cluster",
  label = TRUE,
  repel = TRUE,
  label.size = 40,
  pt.size = 1.0,
  cols = oligo_anno_colors,
  raster = FALSE
) +
  ggtitle(NULL) +
  labs(color = "Cluster") +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 70),
    legend.title = element_text(size = 70),
    plot.title = element_blank()
  )

p_anno_umap <- add_umap_arrows(oligo_data, p_anno_umap, reduction_name)

ggsave(
  file.path(output_dir, "plots/UMAP_oligo_clusters.png"),
  p_anno_umap,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)

# 额外保存一份更明确的文件名
ggsave(
  file.path(output_dir, "plots/UMAP_oligo_annotated_clusters.png"),
  p_anno_umap,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)

# ============================================================================
# 2. UMAP by Sample
# ============================================================================

if ("sample" %in% colnames(oligo_data@meta.data)) {
  samples <- sort(unique(as.character(oligo_data$sample)))
  sample_colors <- setNames(
    colorRampPalette(subcluster_palette_olig)(length(samples)),
    samples
  )

  p_sample_umap <- DimPlot(
    oligo_data,
    reduction = reduction_name,
    group.by = "sample",
    label = FALSE,
    pt.size = 1.0,
    cols = sample_colors,
    raster = FALSE
  ) +
    ggtitle(NULL) +
    labs(color = "Sample") +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 70),
      legend.title = element_text(size = 70),
      plot.title = element_blank()
    )

  p_sample_umap <- add_umap_arrows(oligo_data, p_sample_umap, reduction_name)

  ggsave(
    file.path(output_dir, "plots/UMAP_oligo_by_sample.png"),
    p_sample_umap,
    width = 10,
    height = 8,
    dpi = 300,
    bg = "white"
  )
}

# ============================================================================
# 3. UMAP by Genotype
# ============================================================================

if ("genotype" %in% colnames(oligo_data@meta.data)) {
  genotypes <- sort(unique(as.character(oligo_data$genotype)))
  genotype_colors <- setNames(
    colorRampPalette(subcluster_palette_olig)(length(genotypes)),
    genotypes
  )

  p_genotype_umap <- DimPlot(
    oligo_data,
    reduction = reduction_name,
    group.by = "genotype",
    label = FALSE,
    pt.size = 1.0,
    cols = genotype_colors,
    raster = FALSE
  ) +
    ggtitle(NULL) +
    labs(color = "Genotype") +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 70),
      legend.title = element_text(size = 70),
      plot.title = element_blank()
    )

  p_genotype_umap <- add_umap_arrows(oligo_data, p_genotype_umap, reduction_name)

  ggsave(
    file.path(output_dir, "plots/UMAP_oligo_by_genotype.png"),
    p_genotype_umap,
    width = 10,
    height = 8,
    dpi = 300,
    bg = "white"
  )
}

# ============================================================================
# 4. FeaturePlot
# ============================================================================

key_markers <- c("Mbp", "Plp1", "Mag", "Pdgfra", "Cacng4", "Sox10", "Olig1", "Olig2")
key_markers <- key_markers[key_markers %in% rownames(oligo_data)]

if (length(key_markers) > 0) {
  fp <- FeaturePlot(
    oligo_data,
    features = key_markers,
    reduction = reduction_name,
    ncol = 4,
    pt.size = 0.5,
    order = TRUE,
    cols = c("lightgrey", "#D19246")
  ) +
    plot_annotation(
      title = "Key Marker Expression",
      theme = theme(plot.title = element_text(size = 50, face = "bold", hjust = 0.5))
    )

  ggsave(
    file.path(output_dir, "plots/oligo_key_markers_featureplot.png"),
    fp,
    width = 16,
    height = 4 * ceiling(length(key_markers) / 4),
    dpi = 300,
    bg = "white"
  )
}

# ============================================================================
# 5. DotPlot：使用标注后的 cluster
# ============================================================================

dot_genes <- unique(unlist(list(
  "Mature" = c("Mbp", "Plp1", "Mag"),
  "Immature/OPC" = c("Pdgfra", "Cacng4", "Sox10")
)))

dot_genes <- dot_genes[dot_genes %in% rownames(oligo_data)]

if (length(dot_genes) >= 3) {
  dot <- DotPlot(
    oligo_data,
    features = dot_genes,
    group.by = "oligo_annotated_cluster",
    cols = c("lightgrey", "#D19246"),
    dot.scale = 8,
    scale = TRUE
  ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 50),
      axis.text.y = element_text(size = 50),
      axis.title = element_blank(),
      plot.title = element_blank()
    ) +
    labs(x = NULL, y = NULL)

  ggsave(
    file.path(output_dir, "plots/oligo_markers_dotplot_by_annotated_cluster.png"),
    dot,
    width = max(12, length(dot_genes) * 0.8),
    height = 8,
    dpi = 300,
    bg = "white"
  )
}

# ============================================================================
# 6. 标注后 cluster 统计表
# ============================================================================

cluster_stats <- oligo_data@meta.data %>%
  group_by(oligo_annotated_cluster) %>%
  summarise(
    cell_count = n(),
    percentage = n() / ncol(oligo_data) * 100,
    .groups = "drop"
  ) %>%
  arrange(oligo_annotated_cluster)

write.csv(
  cluster_stats,
  file.path(output_dir, "tables/oligo_annotated_cluster_statistics.csv"),
  row.names = FALSE
)

# 同时保存原始 cluster 统计
original_cluster_stats <- oligo_data@meta.data %>%
  group_by(oligo_clusters) %>%
  summarise(
    cell_count = n(),
    percentage = n() / ncol(oligo_data) * 100,
    .groups = "drop"
  ) %>%
  arrange(suppressWarnings(as.numeric(oligo_clusters)))

write.csv(
  original_cluster_stats,
  file.path(output_dir, "tables/oligo_original_cluster_statistics.csv"),
  row.names = FALSE
)

# ============================================================================
# 7. 按 sample 的标注后 cluster 分布图
# 文件名按你要求保存：
# oligo_cluster_distribution_by_sample.png
# ============================================================================

if ("sample" %in% colnames(oligo_data@meta.data)) {
  sample_cluster_dist <- oligo_data@meta.data %>%
    group_by(sample, oligo_annotated_cluster) %>%
    summarise(cell_count = n(), .groups = "drop") %>%
    group_by(sample) %>%
    mutate(
      total_cells = sum(cell_count),
      percentage = cell_count / total_cells * 100
    ) %>%
    ungroup()

  write.csv(
    sample_cluster_dist,
    file.path(output_dir, "tables/oligo_cluster_distribution_by_sample.csv"),
    row.names = FALSE
  )

  p_sample_dist <- ggplot(
    sample_cluster_dist,
    aes(x = sample, y = percentage, fill = oligo_annotated_cluster)
  ) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = oligo_anno_colors, drop = FALSE) +
    theme_minimal(base_family = "Arial") +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 90, color = "black"),
      axis.text.y = element_text(size = 90, color = "black"),
      legend.text = element_text(size = 90),
      legend.title = element_text(size = 90),
      plot.title = element_blank(),
      axis.title = element_blank()
    ) +
    labs(title = NULL, x = NULL, y = NULL, fill = "Cluster")

  ggsave(
    file.path(output_dir, "plots/oligo_cluster_distribution_by_sample.png"),
    p_sample_dist,
    width = 12,
    height = 8,
    dpi = 300,
    bg = "white"
  )
}

# ============================================================================
# 8. 按 genotype 的标注后 cluster 分布图
# 文件名按你要求保存：
# oligo_cluster_distribution_by_genotype.png
# ============================================================================

if ("genotype" %in% colnames(oligo_data@meta.data)) {
  genotype_cluster_dist <- oligo_data@meta.data %>%
    group_by(genotype, oligo_annotated_cluster) %>%
    summarise(cell_count = n(), .groups = "drop") %>%
    group_by(genotype) %>%
    mutate(
      total_cells = sum(cell_count),
      percentage = cell_count / total_cells * 100
    ) %>%
    ungroup()

  write.csv(
    genotype_cluster_dist,
    file.path(output_dir, "tables/oligo_cluster_distribution_by_genotype.csv"),
    row.names = FALSE
  )

  p_genotype_dist <- ggplot(
    genotype_cluster_dist,
    aes(x = genotype, y = percentage, fill = oligo_annotated_cluster)
  ) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = oligo_anno_colors, drop = FALSE) +
    theme_minimal(base_family = "Arial") +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 90, color = "black"),
      axis.text.y = element_text(size = 90, color = "black"),
      legend.text = element_text(size = 90),
      legend.title = element_text(size = 90),
      plot.title = element_blank(),
      axis.title = element_blank()
    ) +
    labs(title = NULL, x = NULL, y = NULL, fill = "Cluster")

  ggsave(
    file.path(output_dir, "plots/oligo_cluster_distribution_by_genotype.png"),
    p_genotype_dist,
    width = 7,
    height = 8,
    dpi = 300,
    bg = "white"
  )
}

# ============================================================================
# 9. DEG：A4-OLs vs Other
# DEG 前不按基因和细胞过滤
# ============================================================================

cat("\n=== Differential Expression: A4-OLs vs Other ===\n")

if ("RNA" %in% Assays(oligo_data)) {
  DefaultAssay(oligo_data) <- "RNA"
}

# Seurat v5 多 layer 情况下需要 JoinLayers
if (exists("JoinLayers")) {
  oligo_data <- JoinLayers(oligo_data)
}

oligo_data$DEG_group <- ifelse(
  as.character(oligo_data$oligo_annotated_cluster) == "A4-OLs",
  "A4-OLs",
  "Other"
)

oligo_data$DEG_group <- factor(oligo_data$DEG_group, levels = c("A4-OLs", "Other"))

deg_group_table <- table(oligo_data$DEG_group)
cat("DEG group cell numbers:\n")
print(deg_group_table)

if (any(deg_group_table < 3)) {
  warning("One DEG group has fewer than 3 cells. FindMarkers may be unstable.")
}

Idents(oligo_data) <- oligo_data$DEG_group

deg_raw <- FindMarkers(
  oligo_data,
  ident.1 = "A4-OLs",
  ident.2 = "Other",
  features = rownames(oligo_data),
  min.pct = 0,
  logfc.threshold = 0,
  min.diff.pct = -Inf,
  return.thresh = Inf,
  test.use = "wilcox",
  only.pos = FALSE,
  verbose = FALSE
)

deg_raw$gene <- rownames(deg_raw)

# 兼容 Seurat 不同版本的 logFC 列名
fc_col <- intersect(c("avg_log2FC", "avg_logFC"), colnames(deg_raw))[1]

if (is.na(fc_col)) {
  stop("Cannot find avg_log2FC or avg_logFC column in FindMarkers result.")
}

deg <- deg_raw %>%
  mutate(avg_log2FC = .data[[fc_col]]) %>%
  transmute(
    gene = gene,
    pvalue = p_val,
    avg_log2FC = avg_log2FC,
    pct_A4_OLs = pct.1,
    pct_Other = pct.2,
    p_val_adj = p_val_adj
  ) %>%
  arrange(pvalue)

write.csv(
  deg,
  file.path(output_dir, "tables/DEG_A4_OLs_vs_Other_all_genes.csv"),
  row.names = FALSE
)

cat(sprintf("DEG finished. Total genes returned: %d\n", nrow(deg)))

# ============================================================================
# 10. 火山图：A4-OLs vs Other
# 沿用之前火山图格式，同时标注 Egr1
# ============================================================================

sig_cutoff <- 0.05
fc_cutoff <- log2(1.15)

df <- deg %>%
  filter(!is.na(p_val_adj), !is.na(avg_log2FC))

# 为了避免 p_val_adj == 0 导致 Inf，仅用于画图替换；原始表格不改
pos_padj <- df$p_val_adj[df$p_val_adj > 0]

if (length(pos_padj) > 0) {
  min_padj_visible <- min(pos_padj)
} else {
  min_padj_visible <- 1e-300
}

df <- df %>%
  mutate(
    padj_plot = ifelse(p_val_adj <= min_padj_visible, min_padj_visible, p_val_adj),
    significant = p_val_adj <= sig_cutoff & abs(avg_log2FC) >= fc_cutoff,
    direction = case_when(
      significant & avg_log2FC > 0 ~ "Up",
      significant & avg_log2FC < 0 ~ "Down",
      TRUE ~ "Not significant"
    )
  )

sig_up_top3 <- df %>%
  filter(significant, avg_log2FC > 0) %>%
  arrange(p_val_adj, desc(abs(avg_log2FC))) %>%
  slice_head(n = 3)

sig_down_top3 <- df %>%
  filter(significant, avg_log2FC < 0) %>%
  arrange(p_val_adj, desc(abs(avg_log2FC))) %>%
  slice_head(n = 3)

# 额外标注 Egr1，不管是否进入上下调 top3
egr1_label <- df %>%
  filter(tolower(gene) == "egr1") %>%
  slice_head(n = 1)

if (nrow(egr1_label) == 0) {
  warning("Egr1 not found in DEG result.")
}

label_genes <- bind_rows(sig_up_top3, sig_down_top3, egr1_label) %>%
  distinct(gene, .keep_all = TRUE)

write.csv(
  df %>% filter(significant),
  file.path(output_dir, "tables/significant_genes_A4_OLs_vs_Other.csv"),
  row.names = FALSE
)

write.csv(
  label_genes,
  file.path(output_dir, "tables/volcano_labeled_genes_A4_OLs_vs_Other.csv"),
  row.names = FALSE
)

n_up_sig <- sum(df$significant & df$avg_log2FC > 0)
n_down_sig <- sum(df$significant & df$avg_log2FC < 0)

cat(sprintf(
  "Significant DEGs: %d | Up: %d | Down: %d\n",
  sum(df$significant),
  n_up_sig,
  n_down_sig
))

volcano <- ggplot(
  df,
  aes(x = avg_log2FC, y = -log10(padj_plot), color = avg_log2FC)
) +
  geom_point(aes(size = -log10(padj_plot)), alpha = 0.6) +
  scale_color_gradientn(
    colours = c("#3288bd", "#66c2a5", "#ffffbf", "#f46d43", "#9e0142"),
    values = seq(0, 1, 0.2),
    name = "log2 FC"
  ) +
  scale_size(name = "-log10(Padj)") +
  geom_vline(
    xintercept = c(-fc_cutoff, fc_cutoff),
    linetype = 2,
    linewidth = 0.3
  ) +
  geom_hline(
    yintercept = -log10(sig_cutoff),
    linetype = 2,
    linewidth = 0.3
  ) +
  geom_text_repel(
    data = label_genes,
    aes(label = gene),
    size = 22,
    family = "Arial",
    box.padding = 0.5,
    point.padding = 0.3,
    max.overlaps = Inf,
    min.segment.length = 0,
    show.legend = FALSE
  ) +
  labs(
    x = "log2 FC",
    y = "-log10(Padj)",
    title = "A4-OLs vs Other"
  ) +
  theme(
    axis.text = element_text(size = 60, color = "black"),
    axis.title = element_text(size = 70, face = "bold"),
    legend.text = element_text(size = 60),
    legend.title = element_text(size = 70),
    plot.title = element_text(size = 70, face = "bold", hjust = 0.5)
  )

ggsave(
  file.path(output_dir, "plots/volcano_A4_OLs_vs_Other.png"),
  volcano,
  width = 12,
  height = 12,
  dpi = 300,
  bg = "white"
)

# ============================================================================
# 11. GO / KEGG 富集分析
# ============================================================================

enrich_analysis <- function(gene_list, type_label, prefix) {
  gene_list <- unique(gene_list)
  gene_list <- gene_list[!is.na(gene_list)]

  if (length(gene_list) < 10) {
    cat(sprintf("  Skip enrichment: %s, gene number < 10\n", type_label))
    return(NULL)
  }

  entrez_ids <- mapIds(
    org.Mm.eg.db,
    keys = gene_list,
    column = "ENTREZID",
    keytype = "SYMBOL",
    multiVals = "first"
  )

  entrez_ids <- unique(entrez_ids[!is.na(entrez_ids)])

  if (length(entrez_ids) < 10) {
    cat(sprintf("  Skip enrichment: %s, Entrez ID number < 10\n", type_label))
    return(NULL)
  }

  for (ont in c("BP", "CC", "MF")) {
    ego <- enrichGO(
      gene = entrez_ids,
      OrgDb = org.Mm.eg.db,
      ont = ont,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = TRUE,
      minGSSize = 10,
      maxGSSize = 500
    )

    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      write.csv(
        as.data.frame(ego),
        file.path(output_dir, "tables", paste0(prefix, "GO_", ont, "_", type_label, ".csv")),
        row.names = FALSE
      )

      dp <- dotplot(
        ego,
        showCategory = 20,
        title = paste("GO", ont, "-", type_label)
      )

      ggsave(
        file.path(output_dir, "plots", paste0(prefix, "GO_", ont, "_dotplot.png")),
        dp,
        width = 12,
        height = 8,
        dpi = 300,
        bg = "white"
      )
    }
  }

  kk <- enrichKEGG(
    gene = entrez_ids,
    organism = "mmu",
    keyType = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    minGSSize = 10,
    maxGSSize = 500
  )

  if (!is.null(kk) && nrow(as.data.frame(kk)) > 0) {
    write.csv(
      as.data.frame(kk),
      file.path(output_dir, "tables", paste0(prefix, "KEGG_", type_label, ".csv")),
      row.names = FALSE
    )

    dp <- dotplot(
      kk,
      showCategory = 20,
      title = paste("KEGG -", type_label)
    )

    ggsave(
      file.path(output_dir, "plots", paste0(prefix, "KEGG_dotplot.png")),
      dp,
      width = 12,
      height = 8,
      dpi = 300,
      bg = "white"
    )
  }
}

real_sig <- df %>%
  filter(significant)

up_genes <- real_sig %>%
  filter(avg_log2FC > 0) %>%
  pull(gene)

down_genes <- real_sig %>%
  filter(avg_log2FC < 0) %>%
  pull(gene)

cat("\nEnrichment analysis...\n")

enrich_analysis(real_sig$gene, "All_Significant", "all_A4_OLs_vs_Other_")
enrich_analysis(up_genes, "Upregulated", "up_A4_OLs_vs_Other_")
enrich_analysis(down_genes, "Downregulated", "down_A4_OLs_vs_Other_")

# ============================================================================
# 12. 保存对象
# ============================================================================

saveRDS(
  oligo_data,
  file.path(output_dir, "seurat_oligodendrocytes_subset_annotated_A4_OLs.rds"),
  compress = TRUE
)

cat("\nAll results saved to:", output_dir, "\n")

showtext_auto(FALSE)
