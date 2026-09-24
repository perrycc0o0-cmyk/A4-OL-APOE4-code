#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(showtext)
  library(sysfonts)
})

source(file.path("R", "load_config.R"))

if (!file.exists(A4OL_FONT_FILE)) {
  stop("Arial font file not found: ", A4OL_FONT_FILE)
}
sysfonts::font_add("Arial", regular = A4OL_FONT_FILE)
showtext::showtext_auto()
showtext::showtext_opts(dpi = 300)

input_rds <- file.path(
  A4OL_FIGURE_ROOT,
  "olig-c",
  "seurat_oligodendrocytes_subset_annotated_A4_OLs.rds"
)
output_dir <- file.path(
  A4OL_FIGURE_ROOT,
  "SCENIC_A4OL_FORMAL_V9",
  "SCENIC_EGR1_visualization_all",
  "plots"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(input_rds)
DefaultAssay(obj) <- "RNA"

if (!"Egr1" %in% rownames(obj)) {
  stop("Egr1 is absent from the mouse oligodendrocyte object.")
}

reduction_candidates <- c("umap.olig", "umap_harmony", "umap")
reduction_use <- reduction_candidates[
  reduction_candidates %in% names(obj@reductions)
][1]
if (is.na(reduction_use)) {
  stop("No UMAP reduction was found in the mouse oligodendrocyte object.")
}

p <- FeaturePlot(
  obj,
  features = "Egr1",
  reduction = reduction_use,
  order = TRUE,
  pt.size = 0.5,
  cols = c("grey90", "#C43C39")
) +
  labs(title = "Egr1 expression in mouse oligodendrocytes") +
  theme_void(base_family = "Arial") +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )

ggsave(
  file.path(output_dir, "Extended_Data_Fig8a_mouse_Egr1_UMAP.png"),
  p,
  width = 7,
  height = 6,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(output_dir, "Extended_Data_Fig8a_mouse_Egr1_UMAP.pdf"),
  p,
  width = 7,
  height = 6,
  bg = "white"
)
