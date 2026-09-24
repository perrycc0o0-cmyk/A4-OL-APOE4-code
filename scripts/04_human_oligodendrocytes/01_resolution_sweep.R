# ============================================================================
# Oligodendrocytes 亚群分析：多分辨率批量处理
# 分辨率：0.4, 0.6, 0.8, 1.0, 1.2
# 每个分辨率独立聚类、GSEA、可视化，输出至子文件夹
# ============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(fgsea)
library(cowplot)
library(showtext)

source(file.path("R", "load_config.R"))
source(file.path("R", "map_mouse_human_orthologs.R"))

# ------------------------------ 字体设置 ------------------------------------
font_add("Arial", A4OL_FONT_FILE)
showtext_auto()

theme_set(theme_minimal(base_family = "Arial", base_size = 50) +
          theme(panel.grid = element_blank(),
                axis.text = element_text(color = "black", size = 50),
                axis.title = element_text(size = 50),
                plot.title = element_text(size = 50, face = "bold", hjust = 0.5),
                legend.text = element_text(size = 50),
                legend.title = element_text(size = 50)))

# ------------------------------ 路径与参数 --------------------------------
input_rds <- file.path(
  A4OL_PROJECT_ROOT,
  "20260419",
  "snRNAseq-2",
  "integrated_harmony_manual_celltypes.rds"
)
top200_csv <- file.path(
  A4OL_RESULTS_ROOT,
  "Oligodendrocytes_subset",
  "marker",
  "cluster1_vs_others_top200_up_DEGs.csv"
)
base_output_dir <- file.path(A4OL_FIGURE_ROOT, "Mapping", "res")

# 需要处理的分辨率列表
resolutions <- c(0.4, 0.6, 0.8, 1.0, 1.2)
n_pcs <- 20

# UMAP 箭头函数（全局复用）
add_umap_arrows <- function(obj, p, reduction_name) {
  emb <- Embeddings(obj, reduction_name)
  x_range <- range(emb[,1])
  y_range <- range(emb[,2])
  arrow_start_x <- x_range[1] - diff(x_range)*0.05
  arrow_start_y <- y_range[1] - diff(y_range)*0.05
  len_x <- diff(x_range)*0.15
  len_y <- diff(y_range)*0.15
  offset_x <- diff(x_range)*0.03
  offset_y <- diff(y_range)*0.03

  p +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank()) +
    annotate("segment", x = arrow_start_x, xend = arrow_start_x + len_x,
             y = arrow_start_y, yend = arrow_start_y,
             arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
             linewidth = 1.2, color = "black") +
    annotate("segment", x = arrow_start_x, xend = arrow_start_x,
             y = arrow_start_y, yend = arrow_start_y + len_y,
             arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
             linewidth = 1.2, color = "black") +
    annotate("text", x = arrow_start_x + offset_x, y = arrow_start_y - offset_y,
             label = "UMAP-1", hjust = 0, vjust = 1,
             size = 50, family = "Arial") +
    annotate("text", x = arrow_start_x - offset_x,
             y = arrow_start_y + offset_y * 8,
             label = "UMAP-2", hjust = 1, vjust = 0, angle = 90,
             size = 50, family = "Arial") +
    coord_cartesian(clip = "off")
}

# 特征图绘制函数
plot_gene_feature <- function(obj, gene, reduction, save_dir, save_single = TRUE) {
  p <- FeaturePlot(obj, features = gene, reduction = reduction,
                   cols = c("lightgrey", "darkred"), pt.size = 0.8, order = TRUE) +
    ggtitle(gene) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  if (save_single) {
    ggsave(file.path(save_dir, paste0(gene, "_featureplot.png")),
           p, width = 7, height = 6, dpi = 300, bg = "white")
  }
  p_simple <- p +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none")
  return(p_simple)
}

# =========================================================================
# 主循环：对每个分辨率执行完整分析
# =========================================================================
for (target_resolution in resolutions) {
  cat("\n\n========== 开始处理分辨率:", target_resolution, "==========\n")

  # 创建以分辨率命名的子文件夹
  output_dir <- file.path(base_output_dir, paste0("res_", target_resolution))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "marker_plots"), recursive = TRUE, showWarnings = FALSE)

  # ---------- 1. 提取 Oligodendrocytes ----------
  cat("加载数据并提取 Oligodendrocytes...\n")
  scRNA <- readRDS(input_rds)
  DefaultAssay(scRNA) <- "RNA"

  if (!"cell_type_manual" %in% colnames(scRNA@meta.data)) {
    stop("RDS 中缺少 cell_type_manual 列！")
  }

  olig_cells <- subset(scRNA, subset = cell_type_manual == "Oligodendrocytes")
  cat("提取的 Oligodendrocytes 细胞数:", ncol(olig_cells), "\n")
  rm(scRNA); gc()

  # ---------- 2. 重新降维与聚类 ----------
  cat("重新处理：标准化、高变基因、PCA...\n")
  DefaultAssay(olig_cells) <- "RNA"
  olig_cells <- NormalizeData(olig_cells, verbose = FALSE)
  olig_cells <- FindVariableFeatures(olig_cells, nfeatures = 2000, verbose = FALSE)
  olig_cells <- ScaleData(olig_cells, verbose = FALSE)
  olig_cells <- RunPCA(olig_cells, npcs = 50, verbose = FALSE)

  # Harmony 批次校正
  batch_col <- if ("dataset" %in% colnames(olig_cells@meta.data)) "dataset" else "orig.ident"
  cat("使用批次变量:", batch_col, "\n")

  library(harmony)
  harmony_success <- FALSE
  tryCatch({
    olig_cells <- RunHarmony(olig_cells,
                             group.by.vars = batch_col,
                             reduction = "pca",
                             dims = 1:n_pcs,
                             reduction.save = "harmony_olig",
                             verbose = FALSE)
    harmony_success <- TRUE
  }, error = function(e) {
    tryCatch({
      olig_cells <- RunHarmony(olig_cells,
                               group.by.vars = batch_col,
                               reduction.save = "harmony_olig",
                               verbose = FALSE)
      harmony_success <- TRUE
    }, error = function(e2) {})
  })

  reduction_use <- if (harmony_success && "harmony_olig" %in% names(olig_cells@reductions)) "harmony_olig" else "pca"

  # UMAP
  olig_cells <- RunUMAP(olig_cells,
                        reduction = reduction_use,
                        dims = 1:n_pcs,
                        reduction.name = "umap.olig",
                        verbose = FALSE)

  # 聚类
  olig_cells <- FindNeighbors(olig_cells, reduction = reduction_use, dims = 1:n_pcs, verbose = FALSE)
  olig_cells <- FindClusters(olig_cells, resolution = target_resolution, verbose = FALSE)
  olig_cells$olig_clusters <- olig_cells$seurat_clusters
  cat("Olig 亚群聚类数 (res=", target_resolution, "): ", length(unique(olig_cells$olig_clusters)), "\n")
  print(table(olig_cells$olig_clusters))

  # ---------- 3. UMAP 可视化 ----------
  cat("生成 Oligodendrocytes 亚群 UMAP 图...\n")
  umap_use <- "umap.olig"

  p_umap_cluster <- DimPlot(olig_cells, reduction = umap_use, group.by = "olig_clusters",
                            label = TRUE, repel = TRUE, label.size = 45, pt.size = 0.5) +
    ggtitle(paste0("Oligodendrocytes Subclusters (res=", target_resolution, ")")) +
    theme(legend.position = "none")
  p_umap_cluster <- add_umap_arrows(olig_cells, p_umap_cluster, umap_use)
  ggsave(file.path(output_dir, "plots", "UMAP_olig_subclusters.png"),
         p_umap_cluster, width = 12, height = 10, dpi = 300, bg = "white")

  # APOE UMAP
  if ("apoe_group" %in% colnames(olig_cells@meta.data)) {
    apoe_colors <- c("APOE3" = "#1f77b4", "APOE4" = "#ff7f0e")

    p_umap_apoe <- DimPlot(olig_cells, reduction = umap_use, group.by = "apoe_group",
                           cols = apoe_colors, pt.size = 0.5) +
      ggtitle("Oligodendrocytes by APOE Genotype")
    p_umap_apoe <- add_umap_arrows(olig_cells, p_umap_apoe, umap_use)
    ggsave(file.path(output_dir, "plots", "UMAP_olig_APOE.png"),
           p_umap_apoe, width = 12, height = 10, dpi = 300, bg = "white")

    p_umap_apoe_facet <- DimPlot(olig_cells, reduction = umap_use, group.by = "apoe_group",
                                 cols = apoe_colors, pt.size = 0.5) +
      facet_wrap(~ apoe_group) +
      ggtitle("Oligodendrocytes UMAP by APOE Genotype (Facet)")
    p_umap_apoe_facet <- add_umap_arrows(olig_cells, p_umap_apoe_facet, umap_use)
    ggsave(file.path(output_dir, "plots", "UMAP_olig_APOE_facet.png"),
           p_umap_apoe_facet, width = 16, height = 8, dpi = 300, bg = "white")
  }

  # ---------- 4. APOE 占比图 ----------
  if ("apoe_group" %in% colnames(olig_cells@meta.data)) {
    cat("生成 APOE 占比图...\n")
    comp_data <- olig_cells@meta.data %>%
      group_by(olig_clusters, apoe_group) %>%
      summarise(count = n(), .groups = 'drop') %>%
      group_by(olig_clusters) %>%
      mutate(percentage = count / sum(count) * 100)

    write.csv(comp_data, file.path(output_dir, "tables", "olig_cluster_apoe_composition.csv"), row.names = FALSE)

    p_apoe_percent <- ggplot(comp_data, aes(x = olig_clusters, y = percentage, fill = apoe_group)) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_manual(values = c("APOE3" = "#1f77b4", "APOE4" = "#ff7f0e")) +
      labs(x = "Oligodendrocyte Subcluster", y = "Percentage (%)",
           title = "APOE Genotype Composition per Olig Subcluster") +
      theme(legend.title = element_blank())
    ggsave(file.path(output_dir, "plots", "olig_cluster_apoe_percentage.png"),
           p_apoe_percent, width = 10, height = 8, dpi = 300, bg = "white")

    p_apoe_count <- ggplot(comp_data, aes(x = olig_clusters, y = count, fill = apoe_group)) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_manual(values = c("APOE3" = "#1f77b4", "APOE4" = "#ff7f0e")) +
      labs(x = "Oligodendrocyte Subcluster", y = "Number of cells",
           title = "APOE Genotype Count per Olig Subcluster") +
      theme(legend.title = element_blank())
    ggsave(file.path(output_dir, "plots", "olig_cluster_apoe_count.png"),
           p_apoe_count, width = 10, height = 8, dpi = 300, bg = "white")
  }

  # ---------- 5. 加载外部基因集 ----------
  cat("加载外部 Top200 上调基因...\n")
  gene_df <- read.csv(top200_csv)
  if (!"gene" %in% colnames(gene_df)) stop("CSV 中缺少 'gene' 列！")
  mouse_genes <- unique(gene_df$gene)
  ortholog_result <- map_mouse_to_human_orthologs(
    mouse_genes,
    human_expression_genes = rownames(olig_cells)
  )
  human_genes <- ortholog_result$human_genes
  genes_present <- ortholog_result$expression_matched_genes
  write.csv(
    ortholog_result$audit,
    file.path(output_dir, "tables", "mouse_to_human_ortholog_mapping.csv"),
    row.names = FALSE
  )
  cat("原始小鼠基因数:", length(mouse_genes), "同源映射且在人类对象中存在:", length(genes_present), "\n")
  if (length(genes_present) < 5) next  # 跳过此分辨率

  gene_set_list <- list("Cluster1_Up" = genes_present)
  write.csv(data.frame(gene = genes_present),
            file = file.path(output_dir, "tables", "cluster1_up_genes_human.csv"),
            row.names = FALSE)

  # ---------- 6. GSEA 分析 ----------
  cat("计算每个亚群的排序 logFC 并运行 fgsea...\n")
  olig_cells <- JoinLayers(olig_cells, assay = "RNA")
  clusters <- levels(olig_cells$olig_clusters)
  gsea_res_list <- list()

  for (cl in clusters) {
    cat("  分析簇:", cl, "\n")
    fc <- FoldChange(olig_cells,
                     ident.1 = cl,
                     group.by = "olig_clusters",
                     assay = "RNA",
                     slot = "data")
    ranks <- fc$avg_log2FC
    names(ranks) <- rownames(fc)
    ranks <- sort(ranks, decreasing = TRUE)
    ranks <- ranks[!is.na(ranks)]

    common_genes <- intersect(genes_present, names(ranks))
    if (length(common_genes) < 5) {
      cat("   交集过少 (", length(common_genes), ")，跳过。\n")
      next
    }

    set.seed(123)
    fgsea_out <- fgsea(pathways = gene_set_list,
                       stats = ranks,
                       minSize = 5,
                       maxSize = 500,
                       nPermSimple = 10000,
                       eps = 0)
    fgsea_out$cluster <- cl
    gsea_res_list[[cl]] <- fgsea_out
  }

  gsea_all <- bind_rows(gsea_res_list)
  if (nrow(gsea_all) == 0) {
    cat("没有 cluster 通过交集过滤，跳过 GSEA 图。\n")
  } else {
    if ("leadingEdge" %in% colnames(gsea_all)) {
      gsea_all$leadingEdge <- sapply(gsea_all$leadingEdge, function(x) paste(x, collapse = ";"))
    }
    write.csv(gsea_all,
              file = file.path(output_dir, "tables", "GSEA_cluster1_up_per_cluster.csv"),
              row.names = FALSE)

    # 棒棒糖图
    gsea_plot_data <- gsea_all %>%
      filter(!is.na(NES) & !is.na(pval)) %>%
      mutate(
        neg_log10_p = -log10(pval),
        cluster = factor(cluster, levels = clusters)
      )

    p_lollipop <- ggplot(gsea_plot_data, aes(x = cluster, y = NES)) +
      geom_segment(aes(xend = cluster, yend = 0), color = "grey60", linewidth = 0.8) +
      geom_point(aes(size = neg_log10_p, fill = NES), shape = 21, stroke = 0.5) +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
      scale_size_continuous(range = c(4, 12), name = "-log10(p)") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
      labs(x = paste0("Oligodendrocyte Subcluster (res=", target_resolution, ")"),
           y = "Normalized Enrichment Score (NES)",
           title = "Cluster1 Up-regulated Genes Enrichment",
           fill = "NES") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    ggsave(file.path(output_dir, "plots", "GSEA_lollipop.png"),
           p_lollipop, width = 10, height = 6, dpi = 300, bg = "white")
  }

  # ---------- 7. FeaturePlot：Top50 + 指定基因 ----------
  cat("绘制 Top50 和指定基因的 UMAP 特征图...\n")
  top50_human <- head(genes_present, 50)
  plot_list_top50 <- list()
  for (gene in top50_human) {
    p_simple <- plot_gene_feature(olig_cells, gene, umap_use,
                                  file.path(output_dir, "marker_plots"),
                                  save_single = TRUE)
    plot_list_top50[[gene]] <- p_simple
  }

  if (length(plot_list_top50) > 0) {
    ncol_combine <- 8
    nrow_combine <- ceiling(length(plot_list_top50) / ncol_combine)
    combined_plot <- cowplot::plot_grid(plotlist = plot_list_top50,
                                        ncol = ncol_combine,
                                        align = "hv", axis = "none")
    title_label <- cowplot::ggdraw() +
      cowplot::draw_label("Top 50 Up-regulated Genes (Cluster1) Expression UMAP",
                          fontface = "bold", size = 50)
    final_plot <- cowplot::plot_grid(title_label, combined_plot,
                                     ncol = 1, rel_heights = c(0.05, 1))
    ggsave(file.path(output_dir, "plots", "top50_combined_featureplot.png"),
           final_plot, width = 22, height = 2.8 * nrow_combine + 1,
           dpi = 300, bg = "white", limitsize = FALSE)
  }

  # 指定的 5 个基因
  special_genes <- c("TMSB4X", "CALM1", "NRGN", "YWHAH", "PPIA")
  special_genes <- special_genes[special_genes %in% rownames(olig_cells)]
  for (gene in special_genes) {
    cat("  绘制指定基因:", gene, "\n")
    plot_gene_feature(olig_cells, gene, umap_use,
                      file.path(output_dir, "marker_plots"),
                      save_single = TRUE)
  }

  # ---------- 8. 保存对象 ----------
  saveRDS(olig_cells, file = file.path(output_dir, paste0("oligodendrocytes_res", target_resolution, ".rds")),
          compress = TRUE)
  write.csv(olig_cells@meta.data, file = file.path(output_dir, "tables", "olig_metadata.csv"))

  cat("分辨率", target_resolution, "完成！输出至", output_dir, "\n")
}

cat("\n所有分辨率分析完成！\n")
