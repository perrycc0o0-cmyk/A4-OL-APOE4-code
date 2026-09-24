# ============================================================================
# Oligodendrocytes res=0.4 亚群：Cluster2+12 vs Others 差异分析 + 火山图 + 富集（人类版，严格样式）
# 输出：Figure/Mapping/DEG_cluster2_12
# 修改：DEG 前不进行 gene 筛选和细胞筛选
#
# 火山图修改：
# 1. Top gene 沿用原来的形式：上调 Top3 + 下调 Top3，只标注 gene 名称
# 2. Top gene 名称颜色改为比 EGR1 点颜色更深的深粉色：#B94E59
# 3. EGR1 单独用 #F7A6AC 体现
# 4. EGR1 在火山图中显示到 padj = 1e-5，即 -log10(Padj)=5
# ============================================================================

options(timeout = 180)

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)      # 人类数据库
library(enrichplot)
library(AnnotationDbi)

source(file.path("R", "load_config.R"))

# ---------- 字体设置（与参考代码完全一致）----------
library(showtext)
font_add("Arial", A4OL_FONT_FILE)
showtext_auto()

theme_set(theme_minimal(base_family = "Arial", base_size = 50) +
          theme(panel.grid = element_blank(),
                axis.text = element_text(color = "black"),
                plot.title = element_text(hjust = 0.5, face = "bold", size = 50),
                legend.text = element_text(family = "Arial", size = 50),
                legend.title = element_text(family = "Arial", size = 50)))

# ---------- 路径 ----------
input_rds <- file.path(
  A4OL_FIGURE_ROOT,
  "Mapping",
  "res",
  "res_0.4",
  "oligodendrocytes_res0.4.rds"
)

output_dir <- file.path(A4OL_FIGURE_ROOT, "Mapping", "DEG_cluster2_12")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

# ---------- 火山图突出颜色 ----------
egr1_color <- "#F7A6AC"
top_gene_label_color <- "#B94E59"

# ---------- 加载数据 ----------
cat("Loading oligodendrocytes res=0.4 data...\n")

oligo <- readRDS(input_rds)
DefaultAssay(oligo) <- "RNA"

if (!"olig_clusters" %in% colnames(oligo@meta.data)) {
  stop("RDS 中缺少 olig_clusters 列！")
}

# Seurat v5 多 layer 情况下需要 JoinLayers
if (exists("JoinLayers")) {
  oligo <- JoinLayers(oligo, assay = "RNA")
}

# ---------- 设置分组：cluster2+12 vs 其余 ----------
oligo$group <- ifelse(
  as.character(oligo$olig_clusters) %in% c("2", "12"),
  "Cluster2_12",
  "Others"
)

oligo$group <- factor(oligo$group, levels = c("Cluster2_12", "Others"))

Idents(oligo) <- oligo$group

cat("分组统计:\n")
print(table(oligo$group))

# ---------- 差异分析 ----------
cat("\nDifferential expression: Cluster2_12 vs Others...\n")

deg <- FindMarkers(
  oligo,
  ident.1 = "Cluster2_12",
  ident.2 = "Others",
  features = rownames(oligo),
  min.pct = 0,
  logfc.threshold = 0,
  min.diff.pct = -Inf,
  return.thresh = Inf,
  max.cells.per.ident = Inf,
  min.cells.feature = 0,
  min.cells.group = 0,
  test.use = "wilcox",
  only.pos = FALSE,
  verbose = TRUE
)

deg$gene <- rownames(deg)

# 兼容 Seurat 不同版本的 logFC 列名
fc_col <- intersect(c("avg_log2FC", "avg_logFC"), colnames(deg))[1]

if (is.na(fc_col)) {
  stop("FindMarkers 结果中没有 avg_log2FC 或 avg_logFC 列！")
}

deg$avg_log2FC_use <- deg[[fc_col]]

deg <- deg[, c("gene", "p_val", "avg_log2FC_use", "pct.1", "pct.2", "p_val_adj")]
colnames(deg) <- c("gene", "pvalue", "avg_log2FC", "pct_Cluster2_12", "pct_Other", "p_val_adj")

deg <- deg[order(deg$pvalue), ]

write.csv(
  deg,
  file.path(output_dir, "tables/DEG_cluster2_12_vs_others.csv"),
  row.names = FALSE
)

# ---------- 火山图数据准备 ----------
sig_cutoff <- 0.05
fc_cutoff <- log2(1.15)   # ≈ 0.2

df <- deg[!is.na(deg$p_val_adj), ]

# 极小 padj 替换，避免 Inf
pos_padj <- df$p_val_adj[df$p_val_adj > 0]

if (length(pos_padj) > 0) {
  max_logp <- max(-log10(pos_padj))
  min_padj_visible <- 10^(-max_logp)
} else {
  max_logp <- 5
  min_padj_visible <- 10^(-max_logp)
}

df$p_val_adj[df$p_val_adj <= min_padj_visible] <- min_padj_visible

# 标记显著基因
df$significant <- df$p_val_adj <= sig_cutoff & abs(df$avg_log2FC) >= fc_cutoff

# 计算全部显著基因数
n_up <- sum(df$significant & df$avg_log2FC > 0)
n_down <- sum(df$significant & df$avg_log2FC < 0)
n_total <- n_up + n_down

cat(sprintf("Significant DEGs: %d (Up: %d, Down: %d)\n", n_total, n_up, n_down))

# ---------- 火山图标注基因 ----------
# Top gene 沿用原来的形式：上调 Top3
sig_up <- df %>%
  filter(significant, avg_log2FC > 0) %>%
  arrange(p_val_adj) %>%
  head(3)

# Top gene 沿用原来的形式：下调 Top3
sig_down <- df %>%
  filter(significant, avg_log2FC < 0) %>%
  arrange(p_val_adj) %>%
  head(3)

# EGR1 单独提取
egr1_gene <- df %>%
  filter(toupper(gene) == "EGR1") %>%
  head(1)

if (nrow(egr1_gene) == 0) {
  warning("EGR1 not found in DEG result.")
}

# Top gene 标签：原来的 Top3 上调 + Top3 下调
# 如果 EGR1 正好在 Top3 中，先从 Top gene 标签里去掉，避免重复标注
sig_genes_top <- bind_rows(sig_up, sig_down) %>%
  distinct(gene, .keep_all = TRUE) %>%
  filter(toupper(gene) != "EGR1")

# EGR1 显示位置：只改变火山图显示，不改变 DEG 表中的真实 p_val_adj
if (nrow(egr1_gene) > 0) {
  egr1_gene_plot <- egr1_gene
  egr1_gene_plot$p_val_adj_for_plot <- 1e-5
} else {
  egr1_gene_plot <- egr1_gene
  egr1_gene_plot$p_val_adj_for_plot <- numeric(0)
}

# 保存显著基因列表
write.csv(
  df[df$significant, ],
  file.path(output_dir, "tables/significant_genes_cluster2_12.csv"),
  row.names = FALSE
)

# 保存火山图最终标注基因
sig_genes_label_save <- bind_rows(
  sig_genes_top,
  egr1_gene_plot
) %>%
  distinct(gene, .keep_all = TRUE)

write.csv(
  sig_genes_label_save,
  file.path(output_dir, "tables/volcano_labeled_genes_cluster2_12.csv"),
  row.names = FALSE
)

# ---------- 火山图（严格参考样式）----------
cat("Generating volcano plot...\n")

volcano <- ggplot(
  df,
  aes(x = avg_log2FC, y = -log10(p_val_adj), color = avg_log2FC)
) +
  geom_point(aes(size = -log10(p_val_adj)), alpha = 0.6) +
  scale_color_gradientn(
    colours = c("#3288bd", "#66c2a5", "#ffffbf", "#f46d43", "#9e0142"),
    values = seq(0, 1, 0.2),
    name = "log2 FC"
  ) +
  scale_size(name = "-log10(Padj)") +

  # EGR1：单独用 #F7A6AC 体现，显示到 padj = 1e-5 的水平
  geom_point(
    data = egr1_gene_plot,
    aes(x = avg_log2FC, y = -log10(p_val_adj_for_plot)),
    inherit.aes = FALSE,
    color = egr1_color,
    size = 6.0,
    alpha = 1
  ) +

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

  # Top3 上调 + Top3 下调 gene 名称：比 EGR1 点颜色更深
  geom_text_repel(
    data = sig_genes_top,
    aes(
      x = avg_log2FC,
      y = -log10(p_val_adj),
      label = gene
    ),
    inherit.aes = FALSE,
    color = top_gene_label_color,
    size = 22,
    family = "Arial",
    box.padding = 0.5,
    point.padding = 0.3,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +

  # EGR1 gene 名称：用 #F7A6AC，位置同样显示到 padj = 1e-5
  geom_text_repel(
    data = egr1_gene_plot,
    aes(
      x = avg_log2FC,
      y = -log10(p_val_adj_for_plot),
      label = gene
    ),
    inherit.aes = FALSE,
    color = egr1_color,
    size = 22,
    family = "Arial",
    box.padding = 0.5,
    point.padding = 0.3,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +

  labs(
    x = "log2 FC",
    y = "-log10(Padj)",
    title = "Cluster2+12 vs Other Clusters"
  ) +
  theme(
    axis.text = element_text(size = 60, color = "black"),
    axis.title = element_text(size = 70, face = "bold"),
    legend.text = element_text(size = 60),
    legend.title = element_text(size = 70),
    plot.title = element_text(size = 70, face = "bold", hjust = 0.5)
  )

ggsave(
  file.path(output_dir, "plots/volcano_cluster2_12_vs_others.png"),
  volcano,
  width = 12,
  height = 12,
  dpi = 300,
  bg = "white"
)

# ---------- 富集分析函数（人类数据库）----------
enrich_analysis <- function(gene_list, type_label, prefix) {
  gene_list <- unique(gene_list)
  gene_list <- gene_list[!is.na(gene_list)]

  if (length(gene_list) < 10) {
    cat(sprintf("  Skip enrichment for %s: gene number < 10\n", type_label))
    return()
  }

  entrez_ids <- mapIds(
    org.Hs.eg.db,
    keys = gene_list,
    column = "ENTREZID",
    keytype = "SYMBOL",
    multiVals = "first"
  )

  entrez_ids <- unique(entrez_ids[!is.na(entrez_ids)])

  if (length(entrez_ids) < 10) {
    cat(sprintf("  Insufficient genes for enrichment: %s\n", type_label))
    return()
  }

  for (ont in c("BP", "CC", "MF")) {
    ego <- enrichGO(
      gene = entrez_ids,
      OrgDb = org.Hs.eg.db,
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
    organism = "hsa",
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

# ---------- 执行富集分析 ----------
real_sig <- df[df$significant, ]

up_genes <- real_sig[real_sig$avg_log2FC > 0, "gene"]
down_genes <- real_sig[real_sig$avg_log2FC < 0, "gene"]

cat("\nEnrichment analysis for all significant genes...\n")
enrich_analysis(real_sig$gene, "All_Significant", "all_")

cat("Enrichment for upregulated genes...\n")
enrich_analysis(up_genes, "Upregulated", "up_")

cat("Enrichment for downregulated genes...\n")
enrich_analysis(down_genes, "Downregulated", "down_")

cat("\nAll results saved to", output_dir, "\n")

showtext_auto(FALSE)
