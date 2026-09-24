# ============================================================================
# res = 0.4 Oligodendrocytes 亚群合并标注 + APOE 组成图 + GSEA_lollipop 重绘
#
# 从已有对象开始：
# /Figure/Mapping/res/res_0.4/oligodendrocytes_res0.4.rds
#
# 合并规则：
# cluster 2 + 12            -> A4-OLs
# cluster 0 + 8             -> Oligo1
# cluster 1 + 3 + 5 + 6     -> Oligo2
# cluster 4                 -> Oligo3
# cluster 11                -> Oligo4
# cluster 7 + 10            -> Oligo5
# cluster 9                 -> Oligo6
#
# 配色：柔绿森林
# A4-OLs 固定为 #D19246
#
# 输出：
# 1. UMAP_olig_subclusters_annotated.png
# 2. olig_cluster_apoe_percentage.png
# 3. olig_cluster_apoe_count.png
# 4. GSEA_lollipop.png
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(showtext)
  library(fgsea)
})

source(file.path("R", "load_config.R"))
source(file.path("R", "map_mouse_human_orthologs.R"))

# ------------------------------ 字体设置 ------------------------------------
font_add(
  "Arial",
  A4OL_FONT_FILE
)
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

# ------------------------------ 路径 ----------------------------------------
input_rds <- file.path(
  A4OL_FIGURE_ROOT,
  "Mapping",
  "res",
  "res_0.4",
  "oligodendrocytes_res0.4.rds"
)

top200_csv <- file.path(
  A4OL_RESULTS_ROOT,
  "Oligodendrocytes_subset",
  "marker",
  "cluster1_vs_others_top200_up_DEGs.csv"
)

output_dir <- file.path(A4OL_FIGURE_ROOT, "Mapping", "res", "res_0.4")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

# ---------------------------- Subcluster 配色：柔绿森林 ----------------------
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

# A4-OLs 用柔绿森林里的橙色，保持一致
oligo_anno_colors <- c(
  "Oligo1" = "#B8DBB3",
  "Oligo2" = "#86BC79",
  "Oligo3" = "#71A682",
  "Oligo4" = "#81989B",
  "Oligo5" = "#B5AF8B",
  "Oligo6" = "#7EA4B6",
  "A4-OLs" = "#D19246"
)

oligo_anno_levels <- c(
  "Oligo1",
  "Oligo2",
  "Oligo3",
  "Oligo4",
  "Oligo5",
  "Oligo6",
  "A4-OLs"
)

# ------------------------------ 加载对象 ------------------------------------
cat("Loading existing res=0.4 oligodendrocyte object...\n")

olig_cells <- readRDS(input_rds)

cat("Cells:", ncol(olig_cells), "\n")
cat("Genes:", nrow(olig_cells), "\n")

if ("RNA" %in% Assays(olig_cells)) {
  DefaultAssay(olig_cells) <- "RNA"
}

# ------------------------------ 确定 cluster 列 -----------------------------
if ("olig_clusters" %in% colnames(olig_cells@meta.data)) {
  cluster_col <- "olig_clusters"
} else if ("seurat_clusters" %in% colnames(olig_cells@meta.data)) {
  cluster_col <- "seurat_clusters"
  olig_cells$olig_clusters <- as.character(olig_cells$seurat_clusters)
} else {
  stop("对象中没有 olig_clusters 或 seurat_clusters 列。")
}

olig_cells$olig_clusters <- as.character(olig_cells[[cluster_col]][, 1])

cluster_levels <- unique(olig_cells$olig_clusters)

if (all(grepl("^[0-9]+$", cluster_levels))) {
  cluster_levels <- as.character(sort(as.numeric(cluster_levels)))
} else {
  cluster_levels <- sort(cluster_levels)
}

olig_cells$olig_clusters <- factor(
  olig_cells$olig_clusters,
  levels = cluster_levels
)

cat("\nOriginal res=0.4 clusters:\n")
print(table(olig_cells$olig_clusters))

# ------------------------------ res=0.4 合并标注 ----------------------------
cluster_to_oligo_anno <- list(
  "A4-OLs" = c("2", "12"),
  "Oligo1" = c("0", "8"),
  "Oligo2" = c("1", "3", "5", "6"),
  "Oligo3" = c("4"),
  "Oligo4" = c("11"),
  "Oligo5" = c("7", "10"),
  "Oligo6" = c("9")
)

olig_cells$olig_anno_cluster <- NA_character_

for (anno in names(cluster_to_oligo_anno)) {
  cl_use <- cluster_to_oligo_anno[[anno]]
  olig_cells$olig_anno_cluster[
    as.character(olig_cells$olig_clusters) %in% cl_use
  ] <- anno
}

unmapped_clusters <- unique(
  as.character(olig_cells$olig_clusters)[is.na(olig_cells$olig_anno_cluster)]
)

if (length(unmapped_clusters) > 0) {
  warning(
    "以下 cluster 没有被映射，将标注为 Unknown: ",
    paste(unmapped_clusters, collapse = ", ")
  )

  olig_cells$olig_anno_cluster[is.na(olig_cells$olig_anno_cluster)] <- "Unknown"

  oligo_anno_colors <- c(
    oligo_anno_colors,
    "Unknown" = "#D9D9D9"
  )

  oligo_anno_levels <- c(
    oligo_anno_levels,
    "Unknown"
  )
}

olig_cells$olig_anno_cluster <- factor(
  olig_cells$olig_anno_cluster,
  levels = oligo_anno_levels
)

cat("\nAnnotated res=0.4 clusters:\n")
print(table(olig_cells$olig_anno_cluster))

# ------------------------------ 保存映射表 ----------------------------------
anno_mapping_df <- data.frame(
  original_cluster = unlist(cluster_to_oligo_anno),
  annotated_cluster = rep(
    names(cluster_to_oligo_anno),
    lengths(cluster_to_oligo_anno)
  ),
  stringsAsFactors = FALSE
)

anno_mapping_df <- anno_mapping_df %>%
  arrange(
    factor(annotated_cluster, levels = oligo_anno_levels),
    suppressWarnings(as.numeric(original_cluster))
  )

write.csv(
  anno_mapping_df,
  file.path(output_dir, "tables", "res0.4_olig_cluster_annotation_mapping.csv"),
  row.names = FALSE
)

anno_count_df <- as.data.frame(table(olig_cells$olig_anno_cluster))
colnames(anno_count_df) <- c("annotated_cluster", "cell_count")

anno_count_df <- anno_count_df %>%
  mutate(
    percentage = cell_count / sum(cell_count) * 100
  )

write.csv(
  anno_count_df,
  file.path(output_dir, "tables", "res0.4_olig_annotated_cluster_counts.csv"),
  row.names = FALSE
)

# ============================================================================
# UMAP 箭头函数
# ============================================================================
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

# ------------------------------ 确定 UMAP reduction --------------------------
if ("umap.olig" %in% names(olig_cells@reductions)) {
  umap_use <- "umap.olig"
} else if ("umap.harmony" %in% names(olig_cells@reductions)) {
  umap_use <- "umap.harmony"
} else if ("umap" %in% names(olig_cells@reductions)) {
  umap_use <- "umap"
} else {
  stop("对象中没有 umap.olig / umap.harmony / umap reduction。")
}

cat("\nUsing UMAP reduction:", umap_use, "\n")

# ============================================================================
# 1. 合并标注后的 UMAP
# ============================================================================

cat("Plotting annotated oligodendrocyte UMAP...\n")

p_umap_anno <- DimPlot(
  olig_cells,
  reduction = umap_use,
  group.by = "olig_anno_cluster",
  label = TRUE,
  repel = TRUE,
  label.size = 45,
  pt.size = 0.5,
  cols = oligo_anno_colors,
  raster = FALSE
) +
  ggtitle(NULL) +
  labs(color = "Subcluster") +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 50),
    legend.title = element_text(size = 50),
    plot.title = element_blank()
  )

p_umap_anno <- add_umap_arrows(olig_cells, p_umap_anno, umap_use)

ggsave(
  file.path(output_dir, "plots", "UMAP_olig_subclusters_annotated.png"),
  p_umap_anno,
  width = 12,
  height = 10,
  dpi = 300,
  bg = "white"
)

# ============================================================================
# 2. APOE percentage / count 图
# x 轴：APOE group
# fill：标注后的 olig subcluster
# ============================================================================

if (!"apoe_group" %in% colnames(olig_cells@meta.data)) {
  stop("对象中没有 apoe_group 列，无法绘制 APOE 组成图。")
}

cat("Plotting APOE composition plots using annotated olig subclusters...\n")

comp_data <- olig_cells@meta.data %>%
  group_by(apoe_group, olig_anno_cluster) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(apoe_group) %>%
  mutate(
    percentage = count / sum(count) * 100
  ) %>%
  ungroup()

write.csv(
  comp_data,
  file.path(output_dir, "tables", "olig_cluster_apoe_composition.csv"),
  row.names = FALSE
)

# ------------------------------ Percentage 图 --------------------------------
p_apoe_percent <- ggplot(
  comp_data,
  aes(x = apoe_group, y = percentage, fill = olig_anno_cluster)
) +
  geom_bar(
    stat = "identity",
    position = "stack",
    width = 0.7
  ) +
  scale_fill_manual(
    values = oligo_anno_colors,
    drop = FALSE
  ) +
  labs(
    x = NULL,
    y = "Percentage (%)",
    fill = "Subcluster",
    title = NULL
  ) +
  theme_minimal(base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 70, color = "black"),
    axis.text.y = element_text(size = 70, color = "black"),
    axis.title.y = element_text(size = 70, face = "bold"),
    axis.title.x = element_blank(),
    legend.text = element_text(size = 60),
    legend.title = element_text(size = 65, face = "bold"),
    plot.title = element_blank()
  )

ggsave(
  file.path(output_dir, "plots", "olig_cluster_apoe_percentage.png"),
  p_apoe_percent,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)

# ------------------------------ Count 图 -------------------------------------
p_apoe_count <- ggplot(
  comp_data,
  aes(x = apoe_group, y = count, fill = olig_anno_cluster)
) +
  geom_bar(
    stat = "identity",
    position = "stack",
    width = 0.7
  ) +
  scale_fill_manual(
    values = oligo_anno_colors,
    drop = FALSE
  ) +
  labs(
    x = NULL,
    y = "Number of cells",
    fill = "Subcluster",
    title = NULL
  ) +
  theme_minimal(base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 70, color = "black"),
    axis.text.y = element_text(size = 70, color = "black"),
    axis.title.y = element_text(size = 70, face = "bold"),
    axis.title.x = element_blank(),
    legend.text = element_text(size = 60),
    legend.title = element_text(size = 65, face = "bold"),
    plot.title = element_blank()
  )

ggsave(
  file.path(output_dir, "plots", "olig_cluster_apoe_count.png"),
  p_apoe_count,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)

# ============================================================================
# 3. GSEA_lollipop.png：按合并后的 olig_anno_cluster 重新计算和绘图
# ============================================================================

cat("Recomputing GSEA by annotated oligodendrocyte subcluster...\n")

if (!file.exists(top200_csv)) {
  stop("找不到 top200_csv 文件：", top200_csv)
}

gene_df <- read.csv(top200_csv, stringsAsFactors = FALSE)

if (!"gene" %in% colnames(gene_df)) {
  stop("CSV 中缺少 'gene' 列：", top200_csv)
}

mouse_genes <- unique(gene_df$gene)
ortholog_result <- map_mouse_to_human_orthologs(
  mouse_genes,
  human_expression_genes = rownames(olig_cells)
)
human_genes <- ortholog_result$human_genes
genes_present <- ortholog_result$expression_matched_genes
write.csv(
  ortholog_result$audit,
  file.path(output_dir, "tables", "mouse_to_human_ortholog_mapping_res0.4.csv"),
  row.names = FALSE
)

cat("Top200 原始基因数:", length(rat_genes), "\n")
cat("转换后在对象中存在的基因数:", length(genes_present), "\n")

write.csv(
  data.frame(gene = genes_present),
  file.path(output_dir, "tables", "cluster1_up_genes_human_res0.4_annotated.csv"),
  row.names = FALSE
)

if (length(genes_present) < 5) {
  warning("Top200 基因在对象中匹配少于 5 个，跳过 GSEA_lollipop 重绘。")
} else {
  gene_set_list <- list("Cluster1_Up" = genes_present)

  # Seurat v5 多 layer 情况下需要 JoinLayers
  if (exists("JoinLayers")) {
    olig_cells <- JoinLayers(olig_cells, assay = "RNA")
  }

  clusters_anno <- levels(droplevels(olig_cells$olig_anno_cluster))
  gsea_res_list <- list()

  for (cl in clusters_anno) {
    cat("  GSEA for annotated cluster:", cl, "\n")

    n_cells_cl <- sum(olig_cells$olig_anno_cluster == cl, na.rm = TRUE)

    if (n_cells_cl < 3) {
      cat("    细胞数少于 3，跳过：", cl, "\n")
      next
    }

    fc <- FoldChange(
      olig_cells,
      ident.1 = cl,
      group.by = "olig_anno_cluster",
      assay = "RNA",
      slot = "data"
    )

    fc_col <- intersect(c("avg_log2FC", "avg_logFC"), colnames(fc))[1]

    if (is.na(fc_col)) {
      stop("FoldChange 结果中没有 avg_log2FC 或 avg_logFC 列。")
    }

    ranks <- fc[[fc_col]]
    names(ranks) <- rownames(fc)

    ranks <- ranks[!is.na(ranks)]
    ranks <- sort(ranks, decreasing = TRUE)

    common_genes <- intersect(genes_present, names(ranks))

    if (length(common_genes) < 5) {
      cat("    gene set 交集过少：", length(common_genes), "，跳过。\n")
      next
    }

    set.seed(123)

    fgsea_out <- fgsea(
      pathways = gene_set_list,
      stats = ranks,
      minSize = 5,
      maxSize = 500,
      nPermSimple = 10000,
      eps = 0
    )

    fgsea_out$cluster <- cl
    fgsea_out$n_cells <- n_cells_cl

    gsea_res_list[[cl]] <- fgsea_out
  }

  gsea_all <- bind_rows(gsea_res_list)

  if (nrow(gsea_all) == 0) {
    warning("没有 annotated cluster 通过 GSEA 过滤，未生成 GSEA_lollipop.png。")
  } else {
    if ("leadingEdge" %in% colnames(gsea_all)) {
      gsea_all$leadingEdge <- sapply(
        gsea_all$leadingEdge,
        function(x) paste(x, collapse = ";")
      )
    }

    write.csv(
      gsea_all,
      file.path(output_dir, "tables", "GSEA_cluster1_up_per_annotated_cluster.csv"),
      row.names = FALSE
    )

    gsea_plot_data <- gsea_all %>%
      filter(!is.na(NES), !is.na(pval)) %>%
      mutate(
        neg_log10_p = -log10(pmax(pval, 1e-300)),
        cluster = factor(cluster, levels = oligo_anno_levels)
      ) %>%
      filter(!is.na(cluster))

    p_lollipop <- ggplot(
      gsea_plot_data,
      aes(x = cluster, y = NES)
    ) +
      geom_segment(
        aes(xend = cluster, yend = 0),
        color = "grey60",
        linewidth = 0.8
      ) +
      geom_point(
        aes(size = neg_log10_p, fill = cluster),
        shape = 21,
        color = "black",
        stroke = 0.5
      ) +
      scale_fill_manual(
        values = oligo_anno_colors,
        drop = FALSE,
        name = "Subcluster"
      ) +
      scale_size_continuous(
        range = c(5, 14),
        name = "-log10(p)"
      ) +
      geom_hline(
        yintercept = 0,
        linetype = "dashed",
        color = "grey40",
        linewidth = 0.5
      ) +
      labs(
        x = NULL,
        y = "Normalized Enrichment Score (NES)",
        title = NULL
      ) +
      theme_minimal(base_family = "Arial") +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          size = 65,
          color = "black"
        ),
        axis.text.y = element_text(size = 65, color = "black"),
        axis.title.y = element_text(size = 70, face = "bold"),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 55),
        legend.title = element_text(size = 60, face = "bold"),
        plot.title = element_blank()
      )

    ggsave(
      file.path(output_dir, "plots", "GSEA_lollipop.png"),
      p_lollipop,
      width = 12,
      height = 8,
      dpi = 300,
      bg = "white"
    )
  }
}

# ============================================================================
# 4. 保存带新注释的对象和 metadata
# ============================================================================

saveRDS(
  olig_cells,
  file = file.path(output_dir, "oligodendrocytes_res0.4_annotated.rds"),
  compress = TRUE
)

write.csv(
  olig_cells@meta.data,
  file.path(output_dir, "tables", "olig_metadata_annotated.csv")
)

cat("\nDone. Outputs saved to:\n")
cat(output_dir, "\n")

showtext_auto(FALSE)
