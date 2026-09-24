# Map mouse gene symbols to human orthologs using babelgene.
# The audit table records mapped, unmapped, and expression-matched genes.

map_mouse_to_human_orthologs <- function(
    mouse_genes,
    human_expression_genes = NULL,
    min_support = 3) {
  if (!requireNamespace("babelgene", quietly = TRUE)) {
    stop(
      "Package 'babelgene' is required for mouse-to-human ortholog mapping. ",
      "Install it with install.packages('babelgene')."
    )
  }

  mouse_genes <- unique(as.character(mouse_genes))
  mouse_genes <- mouse_genes[!is.na(mouse_genes) & nzchar(mouse_genes)]

  mapping <- babelgene::orthologs(
    genes = mouse_genes,
    species = "mouse",
    human = FALSE,
    min_support = min_support,
    top = TRUE
  )

  if (nrow(mapping) == 0) {
    stop("No mouse-to-human orthologs were returned by babelgene.")
  }

  audit <- merge(
    data.frame(mouse_symbol_input = mouse_genes, stringsAsFactors = FALSE),
    mapping,
    by.x = "mouse_symbol_input",
    by.y = "symbol",
    all.x = TRUE,
    sort = FALSE
  )

  human_symbol_col <- c("human_symbol", "symbol_human")
  human_symbol_col <- human_symbol_col[human_symbol_col %in% colnames(audit)]
  if (length(human_symbol_col) == 0) {
    stop(
      "Cannot identify the human-symbol column returned by babelgene. Columns: ",
      paste(colnames(audit), collapse = ", ")
    )
  }
  human_symbol_col <- human_symbol_col[1]

  human_genes <- unique(as.character(audit[[human_symbol_col]]))
  human_genes <- human_genes[!is.na(human_genes) & nzchar(human_genes)]

  if (is.null(human_expression_genes)) {
    expression_matched_genes <- human_genes
    audit$present_in_human_expression_object <- NA
  } else {
    expression_matched_genes <- intersect(human_genes, human_expression_genes)
    audit$present_in_human_expression_object <-
      !is.na(audit[[human_symbol_col]]) &
      audit[[human_symbol_col]] %in% human_expression_genes
  }

  list(
    mapping = mapping,
    audit = audit,
    human_genes = human_genes,
    expression_matched_genes = expression_matched_genes
  )
}
