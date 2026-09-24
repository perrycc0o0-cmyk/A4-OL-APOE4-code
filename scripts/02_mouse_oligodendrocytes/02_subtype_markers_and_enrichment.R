#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# =============================================================================
# A4-OLs / Oligo1 / Oligo2 / Oligo3 四个亚群分析
#
# 对四个亚群分别进行：
#   1. 目标亚群 vs 其余全部 oligodendrocytes 的 DEG
#   2. 上调、下调基因分别进行 GO BP / CC / MF 和 KEGG 富集
#   3. 提取该亚群上调基因前 50 名作为 marker genes
#   4. 绘制 Top50 marker genes 的 UMAP FeaturePlot
#
# 比较方式：每个目标亚群分别与其余三个亚群合并比较
#
# 输入对象：
#   seurat_oligodendrocytes_subset_annotated_A4_OLs.rds
#
# 四群定义：
#   cluster 0     -> Oligo1
#   cluster 2     -> Oligo2
#   cluster 3     -> Oligo3
#   cluster 1、4  -> A4-OLs
#
# 修复1：显式使用 dplyr::select，避免 AnnotationDbi::select 冲突
# 修复2：火山图标签层设置 inherit.aes=FALSE，避免继承 plot_class
# 修复3：SeuratObject v5 的 FetchData() 使用 layer='data'，不再使用已废弃的 slot
#
# 物种：Mouse
# GO 数据库：org.Mm.eg.db
# KEGG organism：mmu
# =============================================================================

options(stringsAsFactors = FALSE)
options(timeout = 1800)
set.seed(20260802)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(showtext)
  library(sysfonts)
  library(scales)
})

source(file.path("R", "load_config.R"))

# =============================================================================
# 0. 参数
# =============================================================================

BASE_DIR <- A4OL_PROJECT_ROOT

INPUT_RDS <- file.path(
  BASE_DIR,
  "result-qin/Figure/olig-c/",
  "seurat_oligodendrocytes_subset_annotated_A4_OLs.rds"
)

OUTPUT_DIR <- file.path(
  BASE_DIR,
  "result-qin/Figure/olig-c/",
  "ALL4_OLIG_SUBCLUSTERS_DEG_GO_KEGG_TOP50_UMAP_20260802"
)

TABLE_DIR <- file.path(OUTPUT_DIR, "tables")
PLOT_DIR <- file.path(OUTPUT_DIR, "plots")
ENRICH_DIR <- file.path(OUTPUT_DIR, "enrichment")
TOP50_DIR <- file.path(OUTPUT_DIR, "top50_marker_UMAP")
INDIVIDUAL_DIR <- file.path(TOP50_DIR, "individual_genes")

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(ENRICH_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TOP50_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(INDIVIDUAL_DIR, recursive = TRUE, showWarnings = FALSE)

FONT_FILE <- A4OL_FONT_FILE

if (!file.exists(INPUT_RDS)) {
  stop("Input RDS not found: ", INPUT_RDS)
}

if (!file.exists(FONT_FILE)) {
  stop("Arial font not found: ", FONT_FILE)
}

font_add("Arial", FONT_FILE)
showtext_auto()
showtext_opts(dpi = 300)

PLOT_FONT <- "Arial"

TARGET_GROUPS <- c(
  "A4-OLs",
  "Oligo1",
  "Oligo2",
  "Oligo3"
)

FOUR_GROUP_LEVELS <- c(
  "Oligo1",
  "Oligo2",
  "Oligo3",
  "A4-OLs"
)

FOUR_GROUP_COLORS <- c(
  "Oligo1" = "#B8DBB3",
  "Oligo2" = "#86BC79",
  "Oligo3" = "#71A682",
  "A4-OLs" = "#D19246"
)

PADJ_CUTOFF <- 0.05
LOG2FC_CUTOFF <- 0.25
MIN_PCT <- 0.05
TOP_MARKER_N <- 50
TOP_ENRICH_N <- 20
SAVE_INDIVIDUAL_FEATUREPLOTS <- TRUE

# =============================================================================
# 1. 辅助函数
# =============================================================================

cat0 <- function(...) {
  cat(..., "\n", sep = "")
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

add_four_group_annotation <- function(obj) {
  meta_names <- colnames(obj@meta.data)

  # 只有当现有注释列确实包含完整四群时才直接采用。
  if ("oligo_annotated_cluster" %in% meta_names) {
    annotation_values <- unique(
      as.character(
        obj$oligo_annotated_cluster
      )
    )

    if (all(FOUR_GROUP_LEVELS %in% annotation_values)) {
      obj$oligo_four_group <- factor(
        as.character(
          obj$oligo_annotated_cluster
        ),
        levels = FOUR_GROUP_LEVELS
      )

      return(obj)
    }
  }

  # 否则按照原始 cluster 编号恢复四群。
  if (!"oligo_clusters" %in% meta_names) {
    stop(
      "Cannot define Oligo1/Oligo2/Oligo3/A4-OLs. ",
      "Metadata lacks a complete oligo_annotated_cluster column ",
      "and lacks oligo_clusters."
    )
  }

  cl <- as.character(
    obj$oligo_clusters
  )

  four_group <- dplyr::case_when(
    cl == "0" ~ "Oligo1",
    cl == "2" ~ "Oligo2",
    cl == "3" ~ "Oligo3",
    cl %in% c("1", "4") ~ "A4-OLs",
    TRUE ~ NA_character_
  )

  if (any(is.na(four_group))) {
    bad_clusters <- sort(
      unique(
        cl[is.na(four_group)]
      )
    )

    stop(
      "Unmapped oligo cluster IDs: ",
      paste(bad_clusters, collapse = ", ")
    )
  }

  obj$oligo_four_group <- factor(
    four_group,
    levels = FOUR_GROUP_LEVELS
  )

  obj
}

detect_logfc_col <- function(df) {
  hit <- c(
    "avg_log2FC",
    "avg_logFC"
  )

  hit <- hit[
    hit %in% colnames(df)
  ]

  if (length(hit) == 0) {
    stop(
      "FindMarkers output lacks avg_log2FC/avg_logFC."
    )
  }

  hit[[1]]
}

map_mouse_symbols_to_entrez <- function(genes) {
  genes <- unique(
    genes[
      !is.na(genes) &
        genes != ""
    ]
  )

  if (length(genes) == 0) {
    return(character(0))
  }

  ids <- AnnotationDbi::mapIds(
    org.Mm.eg.db,
    keys = genes,
    column = "ENTREZID",
    keytype = "SYMBOL",
    multiVals = "first"
  )

  ids <- unique(
    as.character(
      ids[
        !is.na(ids)
      ]
    )
  )

  ids
}

run_go_kegg_enrichment <- function(
  gene_symbols,
  universe_symbols,
  group_name,
  direction_name
) {
  out_group_dir <- file.path(
    ENRICH_DIR,
    group_name,
    direction_name
  )

  dir.create(
    out_group_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  gene_ids <- map_mouse_symbols_to_entrez(
    gene_symbols
  )

  universe_ids <- map_mouse_symbols_to_entrez(
    universe_symbols
  )

  mapping_summary <- tibble(
    group = group_name,
    direction = direction_name,
    input_gene_symbols = length(unique(gene_symbols)),
    mapped_entrez_ids = length(gene_ids),
    universe_gene_symbols = length(unique(universe_symbols)),
    universe_entrez_ids = length(universe_ids)
  )

  write_csv(
    mapping_summary,
    file.path(
      out_group_dir,
      "gene_ID_mapping_summary.csv"
    )
  )

  if (length(gene_ids) < 10) {
    cat0(
      "[Enrichment skipped] ",
      group_name,
      " / ",
      direction_name,
      ": mapped genes < 10"
    )

    return(invisible(NULL))
  }

  # ---------------------------------------------------------------------------
  # GO BP / CC / MF
  # ---------------------------------------------------------------------------
  for (ontology in c("BP", "CC", "MF")) {
    ego <- tryCatch(
      enrichGO(
        gene = gene_ids,
        universe = universe_ids,
        OrgDb = org.Mm.eg.db,
        keyType = "ENTREZID",
        ont = ontology,
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        minGSSize = 10,
        maxGSSize = 500,
        readable = TRUE
      ),
      error = function(e) {
        cat0(
          "[GO warning] ",
          group_name,
          " / ",
          direction_name,
          " / ",
          ontology,
          ": ",
          conditionMessage(e)
        )
        NULL
      }
    )

    if (!is.null(ego)) {
      ego_df <- as.data.frame(ego)

      if (nrow(ego_df) > 0) {
        ego_df <- ego_df %>%
          arrange(p.adjust, pvalue)

        write_csv(
          ego_df,
          file.path(
            out_group_dir,
            paste0(
              "GO_",
              ontology,
              "_all.csv"
            )
          )
        )

        ego_sig_df <- ego_df %>%
          filter(
            p.adjust < PADJ_CUTOFF
          )

        write_csv(
          ego_sig_df,
          file.path(
            out_group_dir,
            paste0(
              "GO_",
              ontology,
              "_significant_padj0.05.csv"
            )
          )
        )

        if (nrow(ego_sig_df) > 0) {
          ego_sig <- ego
          ego_sig@result <- ego_sig_df

          p_go <- enrichplot::dotplot(
            ego_sig,
            showCategory = min(
              TOP_ENRICH_N,
              nrow(ego_sig_df)
            ),
            orderBy = "GeneRatio"
          ) +
            ggtitle(
              paste0(
                group_name,
                " ",
                direction_name,
                " GO ",
                ontology
              )
            ) +
            theme_bw(
              base_family = PLOT_FONT,
              base_size = 13
            ) +
            theme(
              plot.title = element_text(
                face = "bold",
                hjust = 0.5,
                size = 15
              ),
              axis.text.y = element_text(
                size = 10,
                color = "black"
              ),
              axis.text.x = element_text(
                size = 11,
                color = "black"
              )
            )

          ggsave(
            file.path(
              out_group_dir,
              paste0(
                "GO_",
                ontology,
                "_dotplot.png"
              )
            ),
            p_go,
            width = 10,
            height = 7,
            dpi = 300,
            bg = "white"
          )

          ggsave(
            file.path(
              out_group_dir,
              paste0(
                "GO_",
                ontology,
                "_dotplot.pdf"
              )
            ),
            p_go,
            width = 10,
            height = 7,
            bg = "white"
          )
        }
      }
    }
  }

  # ---------------------------------------------------------------------------
  # KEGG
  # ---------------------------------------------------------------------------
  ekegg <- tryCatch(
    enrichKEGG(
      gene = gene_ids,
      universe = universe_ids,
      organism = "mmu",
      keyType = "kegg",
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      qvalueCutoff = 1,
      minGSSize = 10,
      maxGSSize = 500
    ),
    error = function(e) {
      cat0(
        "[KEGG warning] ",
        group_name,
        " / ",
        direction_name,
        ": ",
        conditionMessage(e)
      )
      NULL
    }
  )

  if (!is.null(ekegg)) {
    ekegg <- tryCatch(
      setReadable(
        ekegg,
        OrgDb = org.Mm.eg.db,
        keyType = "ENTREZID"
      ),
      error = function(e) ekegg
    )

    kegg_df <- as.data.frame(ekegg)

    if (nrow(kegg_df) > 0) {
      kegg_df <- kegg_df %>%
        arrange(p.adjust, pvalue)

      write_csv(
        kegg_df,
        file.path(
          out_group_dir,
          "KEGG_all.csv"
        )
      )

      kegg_sig_df <- kegg_df %>%
        filter(
          p.adjust < PADJ_CUTOFF
        )

      write_csv(
        kegg_sig_df,
        file.path(
          out_group_dir,
          "KEGG_significant_padj0.05.csv"
        )
      )

      if (nrow(kegg_sig_df) > 0) {
        ekegg_sig <- ekegg
        ekegg_sig@result <- kegg_sig_df

        p_kegg <- enrichplot::dotplot(
          ekegg_sig,
          showCategory = min(
            TOP_ENRICH_N,
            nrow(kegg_sig_df)
          ),
          orderBy = "GeneRatio"
        ) +
          ggtitle(
            paste0(
              group_name,
              " ",
              direction_name,
              " KEGG"
            )
          ) +
          theme_bw(
            base_family = PLOT_FONT,
            base_size = 13
          ) +
          theme(
            plot.title = element_text(
              face = "bold",
              hjust = 0.5,
              size = 15
            ),
            axis.text.y = element_text(
              size = 10,
              color = "black"
            ),
            axis.text.x = element_text(
              size = 11,
              color = "black"
            )
          )

        ggsave(
          file.path(
            out_group_dir,
            "KEGG_dotplot.png"
          ),
          p_kegg,
          width = 10,
          height = 7,
          dpi = 300,
          bg = "white"
        )

        ggsave(
          file.path(
            out_group_dir,
            "KEGG_dotplot.pdf"
          ),
          p_kegg,
          width = 10,
          height = 7,
          bg = "white"
        )
      }
    }
  }

  invisible(NULL)
}

make_featureplot <- function(
  embedding_df,
  expression_vector,
  gene_name
) {
  plot_df <- embedding_df

  plot_df$expression <- as.numeric(
    expression_vector[
      plot_df$cell
    ]
  )

  plot_df$expression[
    is.na(plot_df$expression)
  ] <- 0

  plot_df <- plot_df %>%
    arrange(expression)

  positive_values <- plot_df$expression[
    plot_df$expression > 0
  ]

  if (length(positive_values) > 0) {
    upper_limit <- as.numeric(
      quantile(
        positive_values,
        probs = 0.99,
        na.rm = TRUE
      )
    )
  } else {
    upper_limit <- 1
  }

  if (
    !is.finite(upper_limit) ||
      upper_limit <= 0
  ) {
    upper_limit <- 1
  }

  ggplot(
    plot_df,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      color = expression
    )
  ) +
    geom_point(
      size = 0.22,
      alpha = 0.85
    ) +
    scale_color_gradientn(
      colours = c(
        "grey92",
        "#D9D2E9",
        "#9E8AC8",
        "#5A3E99",
        "#2E004D"
      ),
      limits = c(
        0,
        upper_limit
      ),
      oob = scales::squish
    ) +
    ggtitle(gene_name) +
    theme_void(
      base_family = PLOT_FONT
    ) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 13,
        margin = margin(b = 2)
      ),
      legend.position = "none",
      plot.margin = margin(
        t = 3,
        r = 3,
        b = 3,
        l = 3
      )
    )
}

# =============================================================================
# 2. 读取与检查对象
# =============================================================================

cat0("Loading: ", INPUT_RDS)

oligo <- readRDS(
  INPUT_RDS
)

cat0(
  "Object class: ",
  paste(class(oligo), collapse = ", ")
)

cat0(
  "Cells: ",
  ncol(oligo),
  "; genes: ",
  nrow(oligo)
)

if (!"RNA" %in% names(oligo@assays)) {
  stop("RNA assay is missing.")
}

DefaultAssay(oligo) <- "RNA"

oligo <- join_layers_safe(
  oligo
)

# 重新生成 RNA data layer，确保 FeaturePlot 与 DEG 使用同一归一化表达。
oligo <- NormalizeData(
  oligo,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

oligo <- add_four_group_annotation(
  oligo
)

cat0("\nFour-group distribution:")
print(
  table(
    oligo$oligo_four_group
  )
)

write.csv(
  as.data.frame(
    table(
      oligo$oligo_four_group
    )
  ),
  file.path(
    TABLE_DIR,
    "four_group_cell_counts.csv"
  ),
  row.names = FALSE
)

Idents(oligo) <- "oligo_four_group"

# =============================================================================
# 3. 确定 UMAP reduction
# =============================================================================

reduction_candidates <- c(
  "umap.olig",
  "umap.harmony",
  "umap"
)

reduction_use <- reduction_candidates[
  reduction_candidates %in% names(oligo@reductions)
][1]

if (is.na(reduction_use)) {
  stop(
    "No UMAP reduction found. Available reductions: ",
    paste(
      names(oligo@reductions),
      collapse = ", "
    )
  )
}

cat0(
  "UMAP reduction: ",
  reduction_use
)

embedding_df <- as.data.frame(
  Embeddings(
    oligo,
    reduction = reduction_use
  )
)

colnames(embedding_df)[1:2] <- c(
  "UMAP_1",
  "UMAP_2"
)

embedding_df$cell <- rownames(
  embedding_df
)

# 四群总览 UMAP
overview_df <- embedding_df %>%
  left_join(
    oligo@meta.data %>%
      rownames_to_column("cell") %>%
      dplyr::select(
        cell,
        oligo_four_group
      ),
    by = "cell"
  )

center_df <- overview_df %>%
  group_by(
    oligo_four_group
  ) %>%
  summarise(
    UMAP_1 = median(
      UMAP_1,
      na.rm = TRUE
    ),
    UMAP_2 = median(
      UMAP_2,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

p_overview <- ggplot(
  overview_df,
  aes(
    x = UMAP_1,
    y = UMAP_2,
    color = oligo_four_group
  )
) +
  geom_point(
    size = 0.35,
    alpha = 0.80
  ) +
  geom_text_repel(
    data = center_df,
    aes(
      label = oligo_four_group
    ),
    color = "black",
    size = 5,
    family = PLOT_FONT,
    fontface = "bold",
    segment.color = NA,
    max.overlaps = Inf
  ) +
  scale_color_manual(
    values = FOUR_GROUP_COLORS,
    drop = FALSE
  ) +
  theme_void(
    base_family = PLOT_FONT
  ) +
  theme(
    legend.position = "none"
  )

ggsave(
  file.path(
    PLOT_DIR,
    "four_group_UMAP.png"
  ),
  p_overview,
  width = 6,
  height = 5.5,
  dpi = 300,
  bg = "white"
)

# =============================================================================
# 4. 对 A4-OLs / Oligo1 / Oligo2 / Oligo3 分别进行分析
# =============================================================================

all_marker_summary <- list()

for (target_group in TARGET_GROUPS) {
  cat0(
    "\n============================================================"
  )

  cat0(
    "Analyzing ",
    target_group,
    " vs all other oligodendrocytes"
  )

  cat0(
    "============================================================"
  )

  group_table_dir <- file.path(
    TABLE_DIR,
    target_group
  )

  group_plot_dir <- file.path(
    PLOT_DIR,
    target_group
  )

  group_top50_dir <- file.path(
    TOP50_DIR,
    target_group
  )

  group_individual_dir <- file.path(
    INDIVIDUAL_DIR,
    target_group
  )

  dir.create(
    group_table_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  dir.create(
    group_plot_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  dir.create(
    group_top50_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  dir.create(
    group_individual_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  # ---------------------------------------------------------------------------
  # DEG：目标亚群 vs 其他所有 oligodendrocytes
  # ---------------------------------------------------------------------------
  deg_raw <- FindMarkers(
    object = oligo,
    ident.1 = target_group,
    ident.2 = NULL,
    assay = "RNA",
    slot = "data",
    test.use = "wilcox",
    only.pos = FALSE,
    min.pct = MIN_PCT,
    logfc.threshold = 0,
    return.thresh = Inf,
    verbose = FALSE
  )

  deg_raw <- deg_raw %>%
    rownames_to_column("gene")

  fc_col <- detect_logfc_col(
    deg_raw
  )

  deg_df <- deg_raw %>%
    transmute(
      gene = gene,
      p_value = p_val,
      p_adj_seurat = p_val_adj,
      p_adj_BH = p.adjust(
        p_val,
        method = "BH"
      ),
      avg_log2FC = .data[[fc_col]],
      pct_target = pct.1,
      pct_others = pct.2,
      pct_difference = pct.1 - pct.2
    ) %>%
    mutate(
      direction = case_when(
        p_adj_BH < PADJ_CUTOFF &
          avg_log2FC >= LOG2FC_CUTOFF ~ "Up",

        p_adj_BH < PADJ_CUTOFF &
          avg_log2FC <= -LOG2FC_CUTOFF ~ "Down",

        TRUE ~ "NS"
      )
    ) %>%
    arrange(
      p_adj_BH,
      desc(abs(avg_log2FC))
    )

  write_csv(
    deg_df,
    file.path(
      group_table_dir,
      paste0(
        target_group,
        "_vs_other_oligodendrocytes_DEG_all.csv"
      )
    )
  )

  up_df <- deg_df %>%
    filter(
      direction == "Up"
    ) %>%
    arrange(
      desc(avg_log2FC),
      p_adj_BH
    )

  down_df <- deg_df %>%
    filter(
      direction == "Down"
    ) %>%
    arrange(
      avg_log2FC,
      p_adj_BH
    )

  write_csv(
    up_df,
    file.path(
      group_table_dir,
      paste0(
        target_group,
        "_upregulated_significant.csv"
      )
    )
  )

  write_csv(
    down_df,
    file.path(
      group_table_dir,
      paste0(
        target_group,
        "_downregulated_significant.csv"
      )
    )
  )

  deg_summary <- tibble(
    group = target_group,
    all_tested_genes = nrow(deg_df),
    upregulated_genes = nrow(up_df),
    downregulated_genes = nrow(down_df),
    padj_cutoff = PADJ_CUTOFF,
    log2FC_cutoff = LOG2FC_CUTOFF
  )

  write_csv(
    deg_summary,
    file.path(
      group_table_dir,
      paste0(
        target_group,
        "_DEG_summary.csv"
      )
    )
  )

  print(deg_summary)

  # ---------------------------------------------------------------------------
  # 火山图
  # ---------------------------------------------------------------------------
  volcano_df <- deg_df %>%
    mutate(
      plot_padj = pmax(
        p_adj_BH,
        .Machine$double.xmin
      ),
      neg_log10_padj = -log10(
        plot_padj
      ),
      plot_class = factor(
        direction,
        levels = c(
          "Down",
          "NS",
          "Up"
        )
      )
    )

  label_df <- bind_rows(
    up_df %>%
      slice_head(n = 8),
    down_df %>%
      slice_head(n = 8)
  ) %>%
    mutate(
      plot_padj = pmax(
        p_adj_BH,
        .Machine$double.xmin
      ),
      neg_log10_padj = -log10(
        plot_padj
      )
    )

  p_volcano <- ggplot(
    volcano_df,
    aes(
      x = avg_log2FC,
      y = neg_log10_padj,
      color = plot_class
    )
  ) +
    geom_point(
      size = 1.2,
      alpha = 0.70
    ) +
    geom_vline(
      xintercept = c(
        -LOG2FC_CUTOFF,
        LOG2FC_CUTOFF
      ),
      linetype = "dashed",
      color = "grey50"
    ) +
    geom_hline(
      yintercept = -log10(
        PADJ_CUTOFF
      ),
      linetype = "dashed",
      color = "grey50"
    ) +
    geom_text_repel(
      data = label_df,
      aes(
        x = avg_log2FC,
        y = neg_log10_padj,
        label = gene
      ),
      inherit.aes = FALSE,
      color = "black",
      size = 3.8,
      family = PLOT_FONT,
      max.overlaps = Inf,
      box.padding = 0.40,
      point.padding = 0.25,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = c(
        "Down" = "#4DBBD5",
        "NS" = "grey75",
        "Up" = "#E64B35"
      ),
      drop = FALSE
    ) +
    labs(
      title = paste0(
        target_group,
        " vs other oligodendrocytes"
      ),
      x = "log2 fold change",
      y = "-log10(BH-adjusted P)",
      color = NULL
    ) +
    theme_classic(
      base_family = PLOT_FONT,
      base_size = 13
    ) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 16
      ),
      axis.title = element_text(
        size = 14
      ),
      axis.text = element_text(
        size = 11,
        color = "black"
      ),
      legend.position = "top"
    )

  ggsave(
    file.path(
      group_plot_dir,
      paste0(
        target_group,
        "_vs_others_volcano.png"
      )
    ),
    p_volcano,
    width = 7.5,
    height = 6.5,
    dpi = 300,
    bg = "white"
  )

  ggsave(
    file.path(
      group_plot_dir,
      paste0(
        target_group,
        "_vs_others_volcano.pdf"
      )
    ),
    p_volcano,
    width = 7.5,
    height = 6.5,
    bg = "white"
  )

  # ---------------------------------------------------------------------------
  # GO / KEGG：上调和下调分别富集
  # ---------------------------------------------------------------------------
  universe_symbols <- deg_df$gene

  run_go_kegg_enrichment(
    gene_symbols = up_df$gene,
    universe_symbols = universe_symbols,
    group_name = target_group,
    direction_name = "Up"
  )

  run_go_kegg_enrichment(
    gene_symbols = down_df$gene,
    universe_symbols = universe_symbols,
    group_name = target_group,
    direction_name = "Down"
  )

  # ---------------------------------------------------------------------------
  # Top50 上调 marker genes
  #
  # 优先选择显著上调基因。
  # 若显著上调基因不足 50 个，则用其余正 logFC 基因补足。
  # ---------------------------------------------------------------------------
  significant_marker_df <- up_df %>%
    mutate(
      marker_source = "significant_up"
    )

  fallback_marker_df <- deg_df %>%
    filter(
      avg_log2FC > 0,
      !gene %in% significant_marker_df$gene
    ) %>%
    arrange(
      desc(avg_log2FC),
      p_adj_BH
    ) %>%
    mutate(
      marker_source = "positive_logFC_fallback"
    )

  top_marker_df <- bind_rows(
    significant_marker_df,
    fallback_marker_df
  ) %>%
    distinct(
      gene,
      .keep_all = TRUE
    ) %>%
    slice_head(
      n = TOP_MARKER_N
    ) %>%
    mutate(
      marker_rank = row_number()
    )

  write_csv(
    top_marker_df,
    file.path(
      group_table_dir,
      paste0(
        target_group,
        "_top50_up_marker_genes.csv"
      )
    )
  )

  all_marker_summary[[target_group]] <- top_marker_df %>%
    dplyr::select(
      marker_rank,
      gene,
      avg_log2FC,
      p_adj_BH,
      pct_target,
      pct_others,
      marker_source
    ) %>%
    mutate(
      group = target_group,
      .before = 1
    )

  marker_genes <- top_marker_df$gene

  marker_genes <- marker_genes[
    marker_genes %in% rownames(oligo)
  ]

  if (length(marker_genes) == 0) {
    cat0(
      "[FeaturePlot skipped] ",
      target_group,
      ": no marker genes found in object"
    )

    next
  }

  expression_df <- FetchData(
    object = oligo,
    vars = marker_genes,
    cells = embedding_df$cell,
    layer = "data"
  )

  expression_df <- as.data.frame(
    expression_df
  )

  rownames(expression_df) <- embedding_df$cell

  featureplot_list <- list()

  for (gene_name in marker_genes) {
    expression_vector <- expression_df[[gene_name]]

    names(expression_vector) <- rownames(
      expression_df
    )

    p_gene <- make_featureplot(
      embedding_df = embedding_df,
      expression_vector = expression_vector,
      gene_name = gene_name
    )

    featureplot_list[[gene_name]] <- p_gene

    if (isTRUE(SAVE_INDIVIDUAL_FEATUREPLOTS)) {
      ggsave(
        file.path(
          group_individual_dir,
          paste0(
            sprintf(
              "%02d",
              match(
                gene_name,
                marker_genes
              )
            ),
            "_",
            gene_name,
            "_UMAP.png"
          )
        ),
        p_gene,
        width = 4.5,
        height = 4,
        dpi = 300,
        bg = "white"
      )
    }
  }

  combined_featureplot <- wrap_plots(
    featureplot_list,
    ncol = 5
  ) +
    plot_annotation(
      title = paste0(
        target_group,
        " top ",
        length(marker_genes),
        " upregulated marker genes"
      ),
      theme = theme(
        text = element_text(
          family = PLOT_FONT
        ),
        plot.title = element_text(
          face = "bold",
          hjust = 0.5,
          size = 22,
          margin = margin(
            b = 10
          )
        )
      )
    )

  feature_height <- max(
    8,
    ceiling(
      length(marker_genes) / 5
    ) * 3.2
  )

  ggsave(
    file.path(
      group_top50_dir,
      paste0(
        target_group,
        "_top50_up_marker_genes_UMAP.png"
      )
    ),
    combined_featureplot,
    width = 18,
    height = feature_height,
    dpi = 300,
    bg = "white",
    limitsize = FALSE
  )

  ggsave(
    file.path(
      group_top50_dir,
      paste0(
        target_group,
        "_top50_up_marker_genes_UMAP.pdf"
      )
    ),
    combined_featureplot,
    width = 18,
    height = feature_height,
    bg = "white",
    limitsize = FALSE
  )
}

# =============================================================================
# 5. 汇总 Top50 marker 表
# =============================================================================

marker_summary_df <- bind_rows(
  all_marker_summary
)

write_csv(
  marker_summary_df,
  file.path(
    TABLE_DIR,
    "ALL4_OLIG_SUBCLUSTERS_top50_up_marker_genes_combined.csv"
  )
)

cat0(
  "\n============================================================"
)

cat0("DONE")

cat0(
  "Output directory: ",
  OUTPUT_DIR
)

cat0(
  "Main outputs:"
)

cat0(
  "  DEG tables: ",
  TABLE_DIR
)

cat0(
  "  Volcano plots: ",
  PLOT_DIR
)

cat0(
  "  GO/KEGG enrichment: ",
  ENRICH_DIR
)

cat0(
  "  Top50 marker UMAPs: ",
  TOP50_DIR
)

cat0(
  "============================================================"
)

showtext_auto(FALSE)
