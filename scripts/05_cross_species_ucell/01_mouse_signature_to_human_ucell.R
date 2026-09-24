#!/usr/bin/env Rscript
# ============================================================================
# STEP1: UCell analysis for mouse A4-OL signature in human oligodendrocytes
#
# 目的：
# 1) 使用小鼠 A4-OLs vs Other 上调 top200 genes 作为 A4-OL signature
# 2) 在人脑 oligodendrocytes res=0.4 对象中计算 UCell score
# 3) 人脑 A4-OLs 定义严格按 marker-0.4.R：cluster 2 + 12 -> A4-OLs
# 4) 输出 UMAP、violin、APOE3/APOE4、sample/donor 层面统计表
#
# Output directory: <A4OL_FIGURE_ROOT>/UCell_A4OL
# ============================================================================

options(stringsAsFactors = FALSE)
options(timeout = 300)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(UCell)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
})

source(file.path("R", "load_config.R"))

# ------------------------------ paths ---------------------------------------
BASE_DIR <- A4OL_SERVER_ROOT
SC_ROOT <- file.path(BASE_DIR, "scRNAseq")
FIGURE_ROOT <- A4OL_FIGURE_ROOT
OUT_DIR <- file.path(FIGURE_ROOT, "UCell_A4OL")
PLOT_DIR <- file.path(OUT_DIR, "plots")
TABLE_DIR <- file.path(OUT_DIR, "tables")
RDS_DIR <- file.path(OUT_DIR, "rds")
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_DIR, recursive = TRUE, showWarnings = FALSE)

# 人脑 olig res=0.4 对象，优先用 marker-0.4.R 已经产生的对象
HUMAN_OLIG_RES04_RDS_CANDIDATES <- c(
  file.path(FIGURE_ROOT, "Mapping/res/res_0.4/oligodendrocytes_res0.4.rds")
)

# 人脑全对象，只有 HUMAN_OLIG_RES04_RDS 不存在时才会尝试使用
HUMAN_FULL_RDS_CANDIDATES <- c(
  file.path(SC_ROOT, "qin-snRNAseq/20260419/snRNAseq-2/integrated_harmony_manual_celltypes.rds"),
  file.path(SC_ROOT, "snRNAseq-2/integrated_harmony_manual_celltypes.rds")
)

# 小鼠 A4-OLs vs Other 上调 genes 来源，优先使用你原先结果表
A4_DEG_CSV_CANDIDATES <- c(
  file.path(FIGURE_ROOT, "olig-c/tables/DEG_A4_OLs_vs_Other_all_genes.csv"),
  file.path(FIGURE_ROOT, "olig-c/tables/A4_OLs_vs_Other_DEG_all_genes.csv"),
  file.path(FIGURE_ROOT, "olig/tables/DEG_cluster1_vs_others.csv"),
  file.path(SC_ROOT, "qin-snRNAseq/result-qin/Oligodendrocytes_subset/marker/cluster1_vs_others_top200_up_DEGs.csv")
)

N_SIGNATURE_GENES <- 200
MIN_SIGNATURE_GENES_PRESENT <- 10
ORTHOLOG_MIN_SUPPORT <- 3
HUMAN_A4_CLUSTERS_RES04 <- c("2", "12")

A4_COLOR <- "#D19246"
GENOTYPE_COLORS <- c("APOE3" = "#FAD9D2", "APOE4" = "#D6ECF2")
OLIGO_ANNOT_COLORS <- c(
  "Oligo1" = "#B8DBB3",
  "Oligo2" = "#86BC79",
  "Oligo3" = "#71A682",
  "Oligo4" = "#81989B",
  "Oligo5" = "#B5AF8B",
  "Oligo6" = "#7EA4B6",
  "A4-OLs" = "#D19246",
  "Unknown" = "#D9D9D9"
)

FONT_FILE <- A4OL_FONT_FILE
if (file.exists(FONT_FILE) && requireNamespace("showtext", quietly = TRUE)) {
  suppressPackageStartupMessages(library(showtext))
  font_add("Arial", FONT_FILE)
  showtext_auto()
  BASE_FAMILY <- "Arial"
} else {
  BASE_FAMILY <- "sans"
}

theme_set(
  theme_minimal(base_family = BASE_FAMILY, base_size = 18) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.text = element_text(color = "black"),
      legend.title = element_text(color = "black")
    )
)

# ------------------------------ helpers -------------------------------------
cat0 <- function(...) cat(..., "\n", sep = "")

find_by_name <- function(root, filename, max_n = 20) {
  if (!dir.exists(root)) return(character(0))
  cmd <- sprintf("find %s -name %s 2>/dev/null", shQuote(root), shQuote(filename))
  res <- tryCatch(system(cmd, intern = TRUE), error = function(e) character(0))
  unique(head(res, max_n))
}

pick_existing_file <- function(candidates, label, search_root = NULL, filename = NULL) {
  candidates <- unique(c(candidates, if (!is.null(search_root) && !is.null(filename)) find_by_name(search_root, filename) else character(0)))
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) {
    stop("Cannot find ", label, ". Tried:\n", paste(candidates, collapse = "\n"))
  }
  hit[1]
}

find_first_col <- function(meta, candidates, required = TRUE, label = "column") {
  hit <- candidates[candidates %in% colnames(meta)]
  if (length(hit) > 0) return(hit[1])
  if (required) {
    stop("Cannot find ", label, ". Candidate columns: ", paste(candidates, collapse = ", "))
  }
  NULL
}

infer_column_by_values <- function(meta, candidates, patterns, required = FALSE, label = "column") {
  hit <- find_first_col(meta, candidates, required = FALSE, label = label)
  if (!is.null(hit)) return(hit)
  for (cn in colnames(meta)) {
    vals <- unique(as.character(meta[[cn]]))
    vals <- vals[!is.na(vals)]
    vals <- head(vals, 200)
    if (any(grepl(patterns, vals, ignore.case = TRUE))) return(cn)
  }
  if (required) stop("Cannot infer ", label)
  NULL
}

get_assay_matrix <- function(obj, assay = "RNA", layer_or_slot = "data") {
  DefaultAssay(obj) <- assay
  out <- tryCatch(GetAssayData(obj, assay = assay, layer = layer_or_slot), error = function(e) NULL)
  if (is.null(out)) out <- tryCatch(GetAssayData(obj, assay = assay, slot = layer_or_slot), error = function(e) NULL)
  if (is.null(out)) stop("Cannot extract assay matrix: assay=", assay, ", layer/slot=", layer_or_slot)
  out
}

match_human_genes_to_object <- function(gene_vec, object_genes) {
  gene_vec <- unique(na.omit(as.character(gene_vec)))
  gene_vec <- gene_vec[gene_vec != ""]
  object_genes <- unique(as.character(object_genes))
  intersect(gene_vec, object_genes)
}

map_mouse_to_human_orthologs <- function(
  mouse_genes,
  human_object_genes,
  min_support = ORTHOLOG_MIN_SUPPORT
) {
  if (!requireNamespace("babelgene", quietly = TRUE)) {
    stop(
      "Package 'babelgene' is required for mouse-to-human ortholog mapping. ",
      "Install it on the server with install.packages('babelgene')."
    )
  }

  mouse_genes <- unique(na.omit(as.character(mouse_genes)))
  mouse_genes <- mouse_genes[mouse_genes != ""]
  human_object_genes <- unique(as.character(human_object_genes))

  if (length(mouse_genes) == 0) {
    stop("No mouse genes were supplied for ortholog mapping.")
  }

  orthologs <- babelgene::orthologs(
    genes = mouse_genes,
    species = "mouse",
    human = FALSE,
    min_support = min_support,
    top = TRUE
  )

  required_cols <- c("symbol", "human_symbol", "support", "support_n")
  missing_cols <- setdiff(required_cols, colnames(orthologs))
  if (length(missing_cols) > 0) {
    stop(
      "Unexpected babelgene output; missing columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  ortholog_table <- orthologs %>%
    transmute(
      mouse_gene = as.character(symbol),
      human_gene = as.character(human_symbol),
      mouse_entrez = as.character(entrez),
      mouse_ensembl = as.character(ensembl),
      human_entrez = as.character(human_entrez),
      human_ensembl = as.character(human_ensembl),
      support = as.character(support),
      support_n = as.integer(support_n)
    )

  mapping <- data.frame(
    mouse_rank = seq_along(mouse_genes),
    mouse_gene = mouse_genes,
    stringsAsFactors = FALSE
  ) %>%
    left_join(ortholog_table, by = "mouse_gene") %>%
    arrange(mouse_rank, desc(support_n))

  mapping$in_human_object <-
    !is.na(mapping$human_gene) &
    mapping$human_gene != "" &
    mapping$human_gene %in% human_object_genes

  mapping$selected_for_ucell <- FALSE
  candidate_idx <- which(mapping$in_human_object)
  if (length(candidate_idx) > 0) {
    keep_idx <- candidate_idx[
      !duplicated(mapping$human_gene[candidate_idx])
    ]
    mapping$selected_for_ucell[keep_idx] <- TRUE
  }

  mapping$mapping_status <- dplyr::case_when(
    is.na(mapping$human_gene) | mapping$human_gene == "" ~
      "no_supported_human_ortholog",
    !mapping$in_human_object ~
      "human_ortholog_not_in_expression_object",
    mapping$selected_for_ucell ~
      "used_in_ucell",
    TRUE ~
      "duplicate_human_ortholog"
  )

  mapping
}

normalize_genotype <- function(x) {
  y <- as.character(x)
  y2 <- ifelse(grepl("APOE.?4|E4|4/4|3/4|4", y, ignore.case = TRUE), "APOE4",
               ifelse(grepl("APOE.?3|E3|3/3|3", y, ignore.case = TRUE), "APOE3", NA_character_))
  factor(y2, levels = c("APOE3", "APOE4"))
}

read_a4_signature_genes <- function(csv_candidates, n_top = 200) {
  sig_file <- pick_existing_file(csv_candidates, "A4-OL DEG/signature csv", search_root = SC_ROOT, filename = "cluster1_vs_others_top200_up_DEGs.csv")
  cat0("Using A4 signature source: ", sig_file)
  df <- read.csv(sig_file, check.names = FALSE)
  gene_col <- intersect(c("gene", "Gene", "genes", "Genes", "SYMBOL", "symbol", "gene_name", "external_gene_name"), colnames(df))[1]
  if (is.na(gene_col)) stop("No gene column found in signature file: ", sig_file)
  fc_col <- intersect(c("avg_log2FC", "avg_log2FC_use", "avg_logFC", "log2FC", "logFC"), colnames(df))[1]
  padj_col <- intersect(c("p_val_adj", "padj", "p.adjust", "FDR", "adj.P.Val"), colnames(df))[1]

  df[[gene_col]] <- as.character(df[[gene_col]])
  df <- df[!is.na(df[[gene_col]]) & df[[gene_col]] != "", , drop = FALSE]

  if (!is.na(fc_col)) {
    df <- df[!is.na(df[[fc_col]]) & df[[fc_col]] > 0, , drop = FALSE]
    if (!is.na(padj_col)) {
      df <- df[order(df[[padj_col]], -df[[fc_col]], na.last = TRUE), , drop = FALSE]
    } else {
      df <- df[order(-df[[fc_col]], na.last = TRUE), , drop = FALSE]
    }
  }

  genes <- unique(df[[gene_col]])
  genes <- head(genes, n_top)
  write.csv(data.frame(rank = seq_along(genes), gene = genes, source_file = sig_file),
            file.path(TABLE_DIR, "A4_signature_genes_raw_top200.csv"), row.names = FALSE)
  genes
}

annotate_human_res04_a4 <- function(obj) {
  if (!"olig_clusters" %in% colnames(obj@meta.data)) {
    if ("seurat_clusters" %in% colnames(obj@meta.data)) {
      obj$olig_clusters <- as.character(obj$seurat_clusters)
    } else {
      stop("Human olig object lacks olig_clusters/seurat_clusters.")
    }
  }
  obj$olig_clusters <- as.character(obj$olig_clusters)
  cluster_to_oligo_anno <- list(
    "A4-OLs" = c("2", "12"),
    "Oligo1" = c("0", "8"),
    "Oligo2" = c("1", "3", "5", "6"),
    "Oligo3" = c("4"),
    "Oligo4" = c("11"),
    "Oligo5" = c("7", "10"),
    "Oligo6" = c("9")
  )
  obj$olig_anno_cluster <- NA_character_
  for (anno in names(cluster_to_oligo_anno)) {
    obj$olig_anno_cluster[obj$olig_clusters %in% cluster_to_oligo_anno[[anno]]] <- anno
  }
  obj$olig_anno_cluster[is.na(obj$olig_anno_cluster)] <- "Unknown"
  levels_use <- c("Oligo1", "Oligo2", "Oligo3", "Oligo4", "Oligo5", "Oligo6", "A4-OLs", "Unknown")
  levels_use <- levels_use[levels_use %in% unique(obj$olig_anno_cluster)]
  obj$olig_anno_cluster <- factor(obj$olig_anno_cluster, levels = levels_use)
  obj$human_A4_like <- ifelse(as.character(obj$olig_anno_cluster) == "A4-OLs", "A4-OLs", "Other OLs")
  obj$human_A4_like <- factor(obj$human_A4_like, levels = c("Other OLs", "A4-OLs"))
  obj
}

prepare_human_olig_res04 <- function() {
  olig_file <- pick_existing_file(HUMAN_OLIG_RES04_RDS_CANDIDATES, "human oligodendrocytes res0.4 RDS", search_root = FIGURE_ROOT, filename = "oligodendrocytes_res0.4.rds")
  cat0("Loading human olig res=0.4 object: ", olig_file)
  obj <- readRDS(olig_file)
  if ("RNA" %in% Assays(obj)) DefaultAssay(obj) <- "RNA"
  obj <- annotate_human_res04_a4(obj)
  obj
}

get_reduction_name <- function(obj) {
  candidates <- c("umap.olig", "umap", "UMAP", "harmony_umap", "umap.harmony")
  hit <- candidates[candidates %in% names(obj@reductions)]
  if (length(hit) == 0) return(NULL)
  hit[1]
}

score_ucell <- function(obj, signature_genes, score_name = "A4_OL_UCell") {
  if ("RNA" %in% Assays(obj)) DefaultAssay(obj) <- "RNA"
  expr <- get_assay_matrix(obj, assay = DefaultAssay(obj), layer_or_slot = "data")
  if (ncol(expr) != ncol(obj)) stop("Expression matrix cell number does not match object.")

  # signature_genes 已经由正式的小鼠-人类同源映射产生；这里只进行精确匹配，
  # 不再使用 toupper() 或大小写不敏感匹配作为跨物种转换。
  signature_genes <- match_human_genes_to_object(
    signature_genes,
    rownames(expr)
  )
  if (length(signature_genes) < MIN_SIGNATURE_GENES_PRESENT) {
    stop(
      "Too few ortholog-mapped signature genes are present in the human object: ",
      length(signature_genes),
      ". Check A4_signature_mouse_to_human_ortholog_mapping.csv."
    )
  }
  cat0("UCell signature genes present in human object: ", length(signature_genes))
  score_df <- UCell::ScoreSignatures_UCell(expr, features = list(A4_OL = signature_genes), maxRank = min(1500, nrow(expr) - 1))
  score_df <- as.data.frame(score_df)
  col_use <- grep("A4_OL", colnames(score_df), value = TRUE)[1]
  if (is.na(col_use)) col_use <- colnames(score_df)[1]
  # 用 meta.data 直接写入 numeric 向量，避免 Seurat [[<- 产生嵌套/非标准列导致后续 dplyr 统计报错
  score_vec <- as.numeric(score_df[colnames(obj), col_use])
  names(score_vec) <- colnames(obj)
  obj@meta.data[[score_name]] <- score_vec
  attr(obj, "A4_signature_genes_used") <- signature_genes
  obj
}

# 更稳健的 Wilcoxon 函数：
# 1) 不依赖 dplyr::cur_data() 的内部 tibble 行索引；
# 2) 自动处理 factor/character/numeric/list/1列data.frame；
# 3) 任一组为空或只有一个分组时返回 NA，不让脚本中断。
safe_wilcox <- function(df, value_col, group_col) {
  if (!value_col %in% colnames(df) || !group_col %in% colnames(df)) return(NA_real_)

  value <- df[[value_col]]
  group <- df[[group_col]]

  if (is.data.frame(value)) value <- value[[1]]
  if (is.data.frame(group)) group <- group[[1]]
  if (is.list(value) && !is.atomic(value)) value <- unlist(value, use.names = FALSE)
  if (is.list(group) && !is.atomic(group)) group <- unlist(group, use.names = FALSE)

  value <- suppressWarnings(as.numeric(value))
  group <- as.character(group)

  n_use <- min(length(value), length(group))
  value <- value[seq_len(n_use)]
  group <- group[seq_len(n_use)]

  keep <- !is.na(value) & !is.na(group) & group != ""
  value <- value[keep]
  group <- factor(group[keep])

  if (length(value) < 2) return(NA_real_)
  if (length(levels(droplevels(group))) != 2) return(NA_real_)
  if (any(table(droplevels(group)) < 1)) return(NA_real_)

  tryCatch(
    suppressWarnings(wilcox.test(value ~ group)$p.value),
    error = function(e) NA_real_
  )
}

# 按一个分组变量逐组做 APOE3 vs APOE4 Wilcoxon，避免 summarise(cur_data()) 兼容问题。
grouped_apoe_wilcox <- function(df, cluster_col, value_col, genotype_col, sample_col = NULL) {
  clusters <- unique(as.character(df[[cluster_col]]))
  clusters <- clusters[!is.na(clusters)]

  out <- lapply(clusters, function(cl) {
    sub <- df[as.character(df[[cluster_col]]) == cl, , drop = FALSE]

    if (is.null(sample_col)) {
      n3 <- sum(as.character(sub[[genotype_col]]) == "APOE3", na.rm = TRUE)
      n4 <- sum(as.character(sub[[genotype_col]]) == "APOE4", na.rm = TRUE)
      data.frame(
        olig_anno_cluster = cl,
        n_APOE3 = n3,
        n_APOE4 = n4,
        p_value = safe_wilcox(sub, value_col, genotype_col),
        stringsAsFactors = FALSE
      )
    } else {
      n3 <- length(unique(as.character(sub[[sample_col]])[as.character(sub[[genotype_col]]) == "APOE3"]))
      n4 <- length(unique(as.character(sub[[sample_col]])[as.character(sub[[genotype_col]]) == "APOE4"]))
      data.frame(
        olig_anno_cluster = cl,
        n_samples_APOE3 = n3,
        n_samples_APOE4 = n4,
        p_value = safe_wilcox(sub, value_col, genotype_col),
        stringsAsFactors = FALSE
      )
    }
  })

  out <- dplyr::bind_rows(out)
  out$p_adj_BH <- p.adjust(out$p_value, method = "BH")
  out
}

plot_umap_score <- function(obj, score_col, group_col = "olig_anno_cluster", file_prefix = "human") {
  red <- get_reduction_name(obj)
  if (is.null(red)) {
    warning("No UMAP reduction found, skip UMAP plots.")
    return(invisible(NULL))
  }
  emb <- as.data.frame(Embeddings(obj, red))
  colnames(emb)[1:2] <- c("UMAP_1", "UMAP_2")
  emb$score <- obj@meta.data[[score_col]]
  emb$group <- obj@meta.data[[group_col]]

  p1 <- ggplot(emb, aes(UMAP_1, UMAP_2, color = score)) +
    geom_point(size = 0.25, alpha = 0.85) +
    scale_color_gradientn(colors = c("#D9D9D9", "#FEE8C8", "#F16913", "#7F0000")) +
    labs(color = score_col, title = paste0(file_prefix, " A4-OL UCell score")) +
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
  ggsave(file.path(PLOT_DIR, paste0(file_prefix, "_UMAP_A4_UCell_score.png")), p1, width = 7, height = 6, dpi = 300, bg = "white")
  ggsave(file.path(PLOT_DIR, paste0(file_prefix, "_UMAP_A4_UCell_score.pdf")), p1, width = 7, height = 6, bg = "white")

  p2 <- ggplot(emb, aes(UMAP_1, UMAP_2, color = group)) +
    geom_point(size = 0.25, alpha = 0.85) +
    scale_color_manual(values = OLIGO_ANNOT_COLORS, drop = FALSE) +
    labs(color = group_col, title = paste0(file_prefix, " annotated OL clusters")) +
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
  ggsave(file.path(PLOT_DIR, paste0(file_prefix, "_UMAP_annotated_clusters.png")), p2, width = 7, height = 6, dpi = 300, bg = "white")
  ggsave(file.path(PLOT_DIR, paste0(file_prefix, "_UMAP_annotated_clusters.pdf")), p2, width = 7, height = 6, bg = "white")
}

# ------------------------------ run -----------------------------------------
cat0("========== Step 1: Read mouse A4-OL signature ==========")
a4_genes_raw <- read_a4_signature_genes(A4_DEG_CSV_CANDIDATES, n_top = N_SIGNATURE_GENES)
cat0("Raw A4 signature genes: ", length(a4_genes_raw))

cat0("\n========== Step 2: Human OL UCell scoring ==========")
human_olig <- prepare_human_olig_res04()
cat0("Human OL cells: ", ncol(human_olig))
cat0("Human annotated clusters:")
print(table(human_olig$olig_anno_cluster))

ortholog_mapping <- map_mouse_to_human_orthologs(
  mouse_genes = a4_genes_raw,
  human_object_genes = rownames(human_olig),
  min_support = ORTHOLOG_MIN_SUPPORT
)

write.csv(
  ortholog_mapping,
  file.path(TABLE_DIR, "A4_signature_mouse_to_human_ortholog_mapping.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    mapping_method = "babelgene::orthologs",
    babelgene_version = as.character(utils::packageVersion("babelgene")),
    source_species = "Mus musculus",
    target_species = "Homo sapiens",
    min_support = ORTHOLOG_MIN_SUPPORT,
    top_supported_match_only = TRUE,
    uppercase_fallback = FALSE,
    stringsAsFactors = FALSE
  ),
  file.path(TABLE_DIR, "A4_signature_ortholog_mapping_manifest.csv"),
  row.names = FALSE
)

human_sig_genes <- ortholog_mapping$human_gene[
  ortholog_mapping$selected_for_ucell
]
human_sig_genes <- unique(human_sig_genes)

write.csv(
  ortholog_mapping %>%
    filter(selected_for_ucell) %>%
    select(
      mouse_rank,
      mouse_gene,
      human_gene,
      support_n,
      support
    ),
  file.path(TABLE_DIR, "A4_signature_genes_used_in_human_UCell.csv"),
  row.names = FALSE
)

cat0("Mouse genes with a supported human ortholog: ", sum(!is.na(ortholog_mapping$human_gene)))
cat0("Orthologs present in the human object and used by UCell: ", length(human_sig_genes))

human_olig <- score_ucell(human_olig, human_sig_genes, "A4_OL_UCell")

human_meta_out <- human_olig@meta.data
human_meta_out$cell_barcode <- rownames(human_meta_out)
write.csv(human_meta_out, file.path(TABLE_DIR, "human_olig_res04_metadata_with_A4_UCell.csv"), row.names = FALSE)
saveRDS(human_olig, file.path(RDS_DIR, "human_olig_res04_with_A4_UCell.rds"), compress = TRUE)

plot_umap_score(human_olig, "A4_OL_UCell", "olig_anno_cluster", "human")

plot_df <- human_olig@meta.data %>%
  mutate(
    olig_anno_cluster = as.factor(olig_anno_cluster),
    human_A4_like = as.factor(human_A4_like)
  )

cluster_summary <- plot_df %>%
  group_by(olig_anno_cluster) %>%
  summarise(
    n_cells = n(),
    mean_A4_UCell = mean(A4_OL_UCell, na.rm = TRUE),
    median_A4_UCell = median(A4_OL_UCell, na.rm = TRUE),
    sd_A4_UCell = sd(A4_OL_UCell, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_A4_UCell))
write.csv(cluster_summary, file.path(TABLE_DIR, "human_cluster_A4_UCell_summary.csv"), row.names = FALSE)

p_vln_cluster <- ggplot(plot_df, aes(x = olig_anno_cluster, y = A4_OL_UCell, fill = olig_anno_cluster)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.2) +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.7, linewidth = 0.2) +
  scale_fill_manual(values = OLIGO_ANNOT_COLORS, drop = FALSE) +
  labs(x = NULL, y = "A4-OL UCell score", title = "Human OL clusters: mouse A4-OL signature") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
ggsave(file.path(PLOT_DIR, "human_A4_UCell_by_annotated_cluster_violin.png"), p_vln_cluster, width = 8.5, height = 6, dpi = 300, bg = "white")
ggsave(file.path(PLOT_DIR, "human_A4_UCell_by_annotated_cluster_violin.pdf"), p_vln_cluster, width = 8.5, height = 6, bg = "white")

p_a4_other <- ggplot(plot_df, aes(x = human_A4_like, y = A4_OL_UCell, fill = human_A4_like)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.2) +
  geom_boxplot(width = 0.14, outlier.shape = NA, alpha = 0.75, linewidth = 0.2) +
  scale_fill_manual(values = c("Other OLs" = "#D9D9D9", "A4-OLs" = A4_COLOR), drop = FALSE) +
  labs(x = NULL, y = "A4-OL UCell score", title = "Human A4-like OLs vs Other OLs") +
  theme(legend.position = "none")
ggsave(file.path(PLOT_DIR, "human_A4_UCell_A4OLs_vs_OtherOLs.png"), p_a4_other, width = 5.5, height = 5.5, dpi = 300, bg = "white")
ggsave(file.path(PLOT_DIR, "human_A4_UCell_A4OLs_vs_OtherOLs.pdf"), p_a4_other, width = 5.5, height = 5.5, bg = "white")

cell_wilcox_a4 <- data.frame(
  comparison = "A4-OLs_vs_Other_OLs",
  p_value = safe_wilcox(plot_df, "A4_OL_UCell", "human_A4_like")
)
write.csv(cell_wilcox_a4, file.path(TABLE_DIR, "human_A4_UCell_A4OLs_vs_OtherOLs_cell_level_wilcox.csv"), row.names = FALSE)

# ------------------------------ APOE genotype -------------------------------
human_genotype_col <- infer_column_by_values(
  human_olig@meta.data,
  candidates = c("apoe_group", "genotype", "APOE", "apoe", "APOE_genotype", "apoe_genotype", "Apoe", "APOE_status"),
  patterns = "APOE.?[34]|E[34]|3/3|3/4|4/4",
  required = FALSE,
  label = "human APOE genotype column"
)
human_sample_col <- find_first_col(
  human_olig@meta.data,
  c("sample", "sample_id", "orig.ident", "donor", "donor_id", "individual", "subject", "patient", "Sample", "library"),
  required = FALSE,
  label = "human sample/donor column"
)

cat0("Detected genotype column: ", ifelse(is.null(human_genotype_col), "NONE", human_genotype_col))
cat0("Detected sample column: ", ifelse(is.null(human_sample_col), "NONE", human_sample_col))

if (!is.null(human_genotype_col)) {
  plot_df$genotype_use <- normalize_genotype(human_olig@meta.data[[human_genotype_col]])
  plot_df$cell_barcode <- rownames(plot_df)

  apoe_cell_summary <- plot_df %>%
    filter(!is.na(genotype_use)) %>%
    group_by(olig_anno_cluster, genotype_use) %>%
    summarise(
      n_cells = n(),
      mean_A4_UCell = mean(A4_OL_UCell, na.rm = TRUE),
      median_A4_UCell = median(A4_OL_UCell, na.rm = TRUE),
      .groups = "drop"
    )
  write.csv(apoe_cell_summary, file.path(TABLE_DIR, "human_A4_UCell_by_cluster_and_APOE_cell_level_summary.csv"), row.names = FALSE)

  apoe_cell_wilcox <- grouped_apoe_wilcox(
    df = plot_df %>% filter(!is.na(genotype_use)),
    cluster_col = "olig_anno_cluster",
    value_col = "A4_OL_UCell",
    genotype_col = "genotype_use",
    sample_col = NULL
  )
  write.csv(apoe_cell_wilcox, file.path(TABLE_DIR, "human_A4_UCell_APOE3_vs_APOE4_cell_level_wilcox_by_cluster.csv"), row.names = FALSE)

  p_apoe_all <- ggplot(plot_df %>% filter(!is.na(genotype_use)), aes(x = genotype_use, y = A4_OL_UCell, fill = genotype_use)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.2) +
    geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.7, linewidth = 0.2) +
    facet_wrap(~ olig_anno_cluster, scales = "free_y") +
    scale_fill_manual(values = GENOTYPE_COLORS, drop = FALSE) +
    labs(x = NULL, y = "A4-OL UCell score", title = "Human OL: A4 signature intensity by APOE") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  ggsave(file.path(PLOT_DIR, "human_A4_UCell_APOE3_vs_APOE4_by_cluster.png"), p_apoe_all, width = 12, height = 7, dpi = 300, bg = "white")
  ggsave(file.path(PLOT_DIR, "human_A4_UCell_APOE3_vs_APOE4_by_cluster.pdf"), p_apoe_all, width = 12, height = 7, bg = "white")

  p_apoe_a4 <- ggplot(plot_df %>% filter(!is.na(genotype_use), olig_anno_cluster == "A4-OLs"), aes(x = genotype_use, y = A4_OL_UCell, fill = genotype_use)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.2) +
    geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.75, linewidth = 0.2) +
    scale_fill_manual(values = GENOTYPE_COLORS, drop = FALSE) +
    labs(x = NULL, y = "A4-OL UCell score", title = "Human A4-like OLs: APOE3 vs APOE4") +
    theme(legend.position = "none")
  ggsave(file.path(PLOT_DIR, "human_A4_like_OL_A4_UCell_APOE3_vs_APOE4.png"), p_apoe_a4, width = 5, height = 5.5, dpi = 300, bg = "white")
  ggsave(file.path(PLOT_DIR, "human_A4_like_OL_A4_UCell_APOE3_vs_APOE4.pdf"), p_apoe_a4, width = 5, height = 5.5, bg = "white")

  if (!is.null(human_sample_col)) {
    plot_df$sample_use <- as.character(human_olig@meta.data[[human_sample_col]])

    sample_score <- plot_df %>%
      filter(!is.na(genotype_use), !is.na(sample_use)) %>%
      group_by(sample_use, genotype_use, olig_anno_cluster) %>%
      summarise(
        n_cells = n(),
        mean_A4_UCell = mean(A4_OL_UCell, na.rm = TRUE),
        median_A4_UCell = median(A4_OL_UCell, na.rm = TRUE),
        .groups = "drop"
      )
    write.csv(sample_score, file.path(TABLE_DIR, "human_sample_level_A4_UCell_by_cluster.csv"), row.names = FALSE)

    sample_level_wilcox <- grouped_apoe_wilcox(
      df = sample_score,
      cluster_col = "olig_anno_cluster",
      value_col = "mean_A4_UCell",
      genotype_col = "genotype_use",
      sample_col = "sample_use"
    )
    write.csv(sample_level_wilcox, file.path(TABLE_DIR, "human_A4_UCell_APOE3_vs_APOE4_sample_level_wilcox_by_cluster.csv"), row.names = FALSE)

    p_sample_score <- ggplot(sample_score, aes(x = genotype_use, y = mean_A4_UCell, fill = genotype_use)) +
      geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.75, linewidth = 0.3) +
      geom_jitter(width = 0.12, size = 1.8, alpha = 0.85) +
      facet_wrap(~ olig_anno_cluster, scales = "free_y") +
      scale_fill_manual(values = GENOTYPE_COLORS, drop = FALSE) +
      labs(x = NULL, y = "Mean A4-OL UCell score per sample", title = "Sample-level A4 signature by APOE") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ggsave(file.path(PLOT_DIR, "human_sample_level_A4_UCell_APOE3_vs_APOE4_by_cluster.png"), p_sample_score, width = 12, height = 7, dpi = 300, bg = "white")
    ggsave(file.path(PLOT_DIR, "human_sample_level_A4_UCell_APOE3_vs_APOE4_by_cluster.pdf"), p_sample_score, width = 12, height = 7, bg = "white")

    comp <- plot_df %>%
      filter(!is.na(genotype_use), !is.na(sample_use)) %>%
      count(sample_use, genotype_use, olig_anno_cluster, name = "n_cluster") %>%
      group_by(sample_use, genotype_use) %>%
      mutate(n_total_OL = sum(n_cluster), prop_cluster = n_cluster / n_total_OL) %>%
      ungroup()
    write.csv(comp, file.path(TABLE_DIR, "human_sample_level_OL_cluster_composition.csv"), row.names = FALSE)

    comp_a4 <- comp %>% filter(olig_anno_cluster == "A4-OLs")
    p_comp <- ggplot(comp_a4, aes(x = genotype_use, y = prop_cluster * 100, fill = genotype_use)) +
      geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.75, linewidth = 0.3) +
      geom_jitter(width = 0.12, size = 1.8, alpha = 0.85) +
      scale_fill_manual(values = GENOTYPE_COLORS, drop = FALSE) +
      labs(x = NULL, y = "A4-like OL proportion among OLs (%)", title = "Sample-level human A4-like OL proportion") +
      theme(legend.position = "none")
    ggsave(file.path(PLOT_DIR, "human_sample_level_A4_like_OL_proportion_APOE3_vs_APOE4.png"), p_comp, width = 5.5, height = 5.5, dpi = 300, bg = "white")
    ggsave(file.path(PLOT_DIR, "human_sample_level_A4_like_OL_proportion_APOE3_vs_APOE4.pdf"), p_comp, width = 5.5, height = 5.5, bg = "white")

    comp_wilcox <- data.frame(
      comparison = "A4_like_OL_proportion_APOE3_vs_APOE4_sample_level",
      p_value = safe_wilcox(comp_a4, "prop_cluster", "genotype_use")
    )
    write.csv(comp_wilcox, file.path(TABLE_DIR, "human_sample_level_A4_like_OL_proportion_APOE3_vs_APOE4_wilcox.csv"), row.names = FALSE)
  }
} else {
  cat0("No APOE genotype column detected. Skip APOE3/APOE4 plots.")
}

cat0("\nDone. Outputs saved to: ", OUT_DIR)
if (exists("showtext_auto")) showtext_auto(FALSE)
