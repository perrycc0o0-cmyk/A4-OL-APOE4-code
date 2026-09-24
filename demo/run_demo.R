#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

input_dir <- file.path("demo", "input")
output_dir <- file.path("demo", "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

markers <- read.delim(
  file.path(input_dir, "mouse_a4_markers.tsv"),
  check.names = FALSE
)
orthologs <- read.delim(
  file.path(input_dir, "mouse_human_orthologs.tsv"),
  check.names = FALSE
)
human_features <- readLines(file.path(input_dir, "human_detected_features.txt"))
human_features <- unique(human_features[nzchar(human_features)])

required_marker_columns <- c("mouse_symbol", "avg_log2FC", "p_val_adj")
required_ortholog_columns <- c("mouse_symbol", "human_symbol", "mapping_status")
stopifnot(all(required_marker_columns %in% colnames(markers)))
stopifnot(all(required_ortholog_columns %in% colnames(orthologs)))

positive <- markers[
  !is.na(markers$avg_log2FC) & markers$avg_log2FC > 0,
  required_marker_columns,
  drop = FALSE
]
positive <- positive[order(-positive$avg_log2FC, positive$p_val_adj), , drop = FALSE]
positive$marker_rank <- seq_len(nrow(positive))

audit <- merge(
  positive,
  orthologs,
  by = "mouse_symbol",
  all.x = TRUE,
  sort = FALSE
)
audit <- audit[match(positive$mouse_symbol, audit$mouse_symbol), , drop = FALSE]
audit$present_in_human_object <-
  !is.na(audit$human_symbol) & audit$human_symbol %in% human_features

human_signature <- unique(
  audit$human_symbol[audit$present_in_human_object]
)

write.table(
  audit,
  file.path(output_dir, "ortholog_mapping_audit.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
writeLines(human_signature, file.path(output_dir, "human_signature.txt"))

message("Mapped mouse markers: ", sum(!is.na(audit$human_symbol)))
message("Human signature genes present: ", length(human_signature))
