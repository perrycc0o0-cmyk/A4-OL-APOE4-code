# ============================================================================
# 细胞注释模块 - 最终版（FeaturePlot 深紫渐变/字体80，Dotplot 字体70）
# ============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(showtext)
library(purrr)

source(file.path("R", "load_config.R"))

# ---------------------------- 字体设置 ----------------------------
font_add("Arial", A4OL_FONT_FILE)
showtext_auto()

theme_set(theme_minimal(base_family = "Arial", base_size = 40) +
          theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 40),
                axis.text = element_text(family = "Arial", size = 40),
                axis.title = element_text(family = "Arial", size = 40),
                legend.text = element_text(family = "Arial", size = 40),
                legend.title = element_text(family = "Arial", size = 40)))

# ---------------------------- 路径与参数 ----------------------------
input_rds <- file.path(A4OL_RESULTS_ROOT, "seurat_clustered.rds")
output_dir <- file.path(A4OL_FIGURE_ROOT, "celltype")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

high_resolution <- 1.0

# ---------------------------- 标记基因（单数名称） ----------------------------
cell_type_markers <- list(
  "Excitatory Neuron"    = c("Slc17a7", "Slc30a3", "Rbfox3", "Map2"),
  "Inhibitory Neuron"    = c("Gad1", "Gad2", "Rbfox3", "Map2"),
  "Oligodendrocyte"      = c("Plp1", "Mbp", "Olig1"),
  "OPC"                  = c("Cacng4", "Pdgfra"),
  "Astrocyte"            = c("Gja1", "Gfap", "Aqp4"),
  "Microglia"            = c("P2ry12", "Cx3cr1"),
  "Vascular"             = c("Cldn5")
)

cluster_to_celltype <- list(
  "Excitatory Neuron"    = c(30,18,10,12,14,42,11,38,8,9,28,3,33,21,1,2,43,32,24),
  "Inhibitory Neuron"    = c(19,16,23,0,6,36,35,22,40,15,7,41,5,20,17),
  "Oligodendrocyte"      = c(4,34,31,39,13),
  "OPC"                  = c(27,37),
  "Microglia"            = c(26,29),
  "Astrocyte"            = c(25),
  "Vascular"             = c()
)

# ---------------------------- Celltype 配色：甜蜜马卡龙 ----------------------------
enhanced_cell_type_colors <- c(
  "Excitatory Neuron"    = "#F7A6AC",  # 粉红
  "Inhibitory Neuron"    = "#EEC78A",  # 玫粉
  "Oligodendrocyte"      = "#B8E5FA",  # 浅蓝
  "OPC"                  = "#B3DDCB",  # 青绿
  "Astrocyte"            = "#CBE4B1",  # 浅绿
  "Microglia"            = "#EEE9A2",  # 浅黄
  "Vascular"             = "#F7B2C7",  # 杏黄
  "Unknown"              = "#D9D9D9"   # 灰色
)
ordered_celltypes <- c("Excitatory Neuron", "Inhibitory Neuron", "Oligodendrocyte",
                       "OPC", "Astrocyte", "Microglia")

# ---------------------------- 加载数据与重新聚类 ----------------------------
cat("Loading data...\n")
scRNA <- readRDS(input_rds)
cat(sprintf("Cells: %d, Genes: %d\n", ncol(scRNA), nrow(scRNA)))

if ("harmony" %in% names(scRNA@reductions)) {
  reduction_use <- "harmony"
} else {
  reduction_use <- "pca"
}
dims_use <- 1:15

scRNA <- FindNeighbors(scRNA, reduction = reduction_use, dims = dims_use, verbose = FALSE)
scRNA <- FindClusters(scRNA, resolution = high_resolution, verbose = FALSE, random.seed = 42)
cluster_column <- paste0("RNA_snn_res.", high_resolution)
scRNA$seurat_clusters_highres <- as.character(scRNA[[cluster_column]][,1])

# 检查标记基因
all_marker_genes <- unique(unlist(cell_type_markers))
available_genes <- all_marker_genes[all_marker_genes %in% rownames(scRNA)]
missing_genes <- setdiff(all_marker_genes, available_genes)
if (length(missing_genes) > 0) {
  cat("Missing markers:", paste(missing_genes, collapse = ", "), "\n")
}

# ---------------------------- 手动标注 ----------------------------
scRNA$cell_type <- "Unknown"
for (ct in names(cluster_to_celltype)) {
  clusters <- as.character(cluster_to_celltype[[ct]])
  exist <- clusters[clusters %in% unique(scRNA$seurat_clusters_highres)]
  if (length(exist) > 0) {
    scRNA$cell_type[scRNA$seurat_clusters_highres %in% exist] <- ct
  }
}

existing_celltypes <- intersect(names(enhanced_cell_type_colors), unique(scRNA$cell_type))
scRNA$cell_type <- factor(scRNA$cell_type, levels = existing_celltypes)

cell_type_mapping <- data.frame(
  cluster = unique(scRNA$seurat_clusters_highres),
  cell_type = sapply(unique(scRNA$seurat_clusters_highres), 
                     function(x) unique(scRNA$cell_type[scRNA$seurat_clusters_highres == x])),
  stringsAsFactors = FALSE
)
write.csv(cell_type_mapping, file.path(output_dir, "tables/cell_type_cluster_mapping.csv"), row.names = FALSE)

# ---------------------------- UMAP 箭头函数 ----------------------------
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
    annotate("segment",
             x = arrow_start_x, xend = arrow_start_x + len_x,
             y = arrow_start_y, yend = arrow_start_y,
             arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
             linewidth = 1.2, color = "black") +
    annotate("segment",
             x = arrow_start_x, xend = arrow_start_x,
             y = arrow_start_y, yend = arrow_start_y + len_y,
             arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
             linewidth = 1.2, color = "black") +
    annotate("text", 
             x = arrow_start_x + offset_x, 
             y = arrow_start_y - offset_y,
             label = "UMAP-1", hjust = 0, vjust = 1, 
             size = 40, family = "Arial") +
    annotate("text", 
             x = arrow_start_x - offset_x, 
             y = arrow_start_y + offset_y * 6,   # UMAP-2 上移
             label = "UMAP-2", hjust = 1, vjust = 0, angle = 90, 
             size = 40, family = "Arial") +
    coord_cartesian(clip = "off")
}

# ---------------------------- 确定 UMAP 降维 ----------------------------
if ("umap.harmony" %in% names(scRNA@reductions)) {
  main_reduction <- "umap.harmony"
} else if ("umap.pca" %in% names(scRNA@reductions)) {
  main_reduction <- "umap.pca"
} else {
  scRNA <- RunUMAP(scRNA, reduction = reduction_use, dims = 1:15, reduction.name = "umap")
  main_reduction <- "umap"
}

# ---------------------------- 自定义 UMAP 细胞类型图（单数、左移、无图例/标题） ----------------------------
umap_emb <- as.data.frame(Embeddings(scRNA, main_reduction))
colnames(umap_emb) <- c("UMAP1", "UMAP2")
umap_emb$cell_type <- scRNA$cell_type

centers <- umap_emb %>%
  group_by(cell_type) %>%
  summarise(x = median(UMAP1), y = median(UMAP2), .groups = 'drop')

x_range <- range(umap_emb$UMAP1)
width <- diff(x_range)
shift <- width / 5

centers <- centers %>%
  mutate(x = ifelse(cell_type == "Excitatory Neuron", x - shift, x))

p_celltype <- ggplot(umap_emb, aes(x = UMAP1, y = UMAP2, color = cell_type)) +
  geom_point(size = 0.3, alpha = 0.6) +
  scale_color_manual(values = enhanced_cell_type_colors, drop = FALSE) +
  geom_text(data = centers, aes(x = x, y = y, label = cell_type), 
            color = "black", size = 120 / ggplot2::.pt, family = "Arial", fontface = "bold",
            inherit.aes = FALSE) +
  theme_void() +
  theme(legend.position = "none",
        plot.title = element_blank())

p_celltype <- add_umap_arrows(scRNA, p_celltype, main_reduction)
ggsave(file.path(output_dir, "plots/UMAP_cell_type_annotation.png"),
       p_celltype, width = 16, height = 12, dpi = 300, bg = "white")

# ---------------------------- 高分辨率聚类 UMAP（保持原样） ----------------------------
n_clusters <- length(unique(scRNA$seurat_clusters_highres))
cluster_colors <- setNames(
  colorRampPalette(c("#1E3A8A","#DC2626","#059669","#7C2D12","#6D28D9",
                     "#0F766E","#BE185D","#374151","#D97706","#0369A1",
                     "#4F46E5","#A16207"))(n_clusters),
  sort(unique(scRNA$seurat_clusters_highres))
)

p_cluster <- DimPlot(scRNA, reduction = main_reduction, 
                     group.by = "seurat_clusters_highres",
                     label = TRUE, repel = TRUE, label.size = 40, pt.size = 0.3,
                     cols = cluster_colors, raster = FALSE) +
  ggtitle(paste0("High-resolution Clusters (Res=", high_resolution, ")")) +
  theme(legend.position = "right",
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 40)) +
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
p_cluster <- add_umap_arrows(scRNA, p_cluster, main_reduction)

ggsave(file.path(output_dir, "plots", paste0("UMAP_clusters_res", high_resolution, ".png")), 
       p_cluster, width = 16, height = 12, dpi = 300, bg = "white")

# ---------------------------- Dotplot（字体70，无标题） ----------------------------
ordered_genes <- c("Slc17a7", "Slc30a3", "Gad1", "Gad2", "Plp1", "Mbp",
                   "Cacng4", "Pdgfra", "Gja1", "Gfap", "P2ry12", "Cx3cr1")
ordered_genes <- ordered_genes[ordered_genes %in% rownames(scRNA)]

scRNA$cell_type <- factor(scRNA$cell_type, levels = ordered_celltypes)

p_dot <- DotPlot(scRNA, features = ordered_genes, group.by = "cell_type",
                 cols = c("lightgrey", "#DC2626"),
                 dot.scale = 14, scale = TRUE, scale.by = "radius", dot.min = 0.05) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 80, face = "bold"),
    axis.text.y = element_text(size = 80, face = "bold"),
    axis.title = element_blank(),
    legend.text = element_text(size = 70),
    legend.title = element_text(size = 70, face = "bold"),
    plot.title = element_blank()
  ) +
  scale_color_gradient2(low = "#1E3A8A", mid = "white", high = "#DC2626", midpoint = 0.5)

ggsave(file.path(output_dir, "plots/key_markers_dotplot_by_celltype.png"),
       p_dot, width = max(18, length(ordered_genes) * 0.45), height = 8, 
       dpi = 300, bg = "white")

# ---------------------------- FeaturePlot（深紫渐变，所有表达细胞染色，字体80） ----------------------------
representative_markers <- c("Slc17a7","Slc30a3","Gad1","Gad2","Plp1","Mbp",
                            "Cacng4","Pdgfra","Gja1","Gfap","P2ry12","Cx3cr1")
representative_markers <- representative_markers[representative_markers %in% rownames(scRNA)]

if (length(representative_markers) > 0) {
  umap_coords <- as.data.frame(Embeddings(scRNA, main_reduction))
  colnames(umap_coords) <- c("UMAP1", "UMAP2")
  umap_coords$cell <- rownames(umap_coords)
  
  expr_data <- FetchData(scRNA, vars = representative_markers)
  plot_data <- cbind(umap_coords, expr_data)
  
  plot_list <- list()
  for (gene in representative_markers) {
    # 直接按表达值着色，使用深紫色渐变
    p <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color = .data[[gene]])) +
      geom_point(size = 0.3) +
      scale_color_gradient(low = "lightgrey", high = "#2E004D", name = NULL) +  # 极深紫色
      ggtitle(gene) +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, size = 80, face = "bold", family = "Arial"),
            legend.position = "none")
    plot_list[[gene]] <- p
  }
  
  combined_feat <- wrap_plots(plot_list, ncol = 4) +
    plot_annotation(title = NULL)
  
  ggsave(file.path(output_dir, "plots/key_marker_genes_feature_plot.png"),
         combined_feat, width = 24, 
         height = 5 * ceiling(length(representative_markers)/4), 
         dpi = 300, bg = "white")
}

# ---------------------------- 细胞比例表 ----------------------------
if ("sample" %in% colnames(scRNA@meta.data) && "genotype" %in% colnames(scRNA@meta.data)) {
  cell_type_composition <- scRNA@meta.data %>%
    group_by(sample, genotype, cell_type) %>%
    summarise(cell_count = n(), .groups = 'drop') %>%
    group_by(sample) %>%
    mutate(percentage = cell_count / sum(cell_count) * 100)
  write.csv(cell_type_composition, 
            file.path(output_dir, "tables/cell_type_composition.csv"), row.names = FALSE)
}

saveRDS(scRNA, file.path(output_dir, "seurat_with_cell_type_annotation.rds"), compress = TRUE)

cat("\nAll outputs saved to:", output_dir, "\n")
showtext_auto(FALSE)
