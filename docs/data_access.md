# Data access

The mouse single-nucleus RNA-sequencing and spatial transcriptomics data
generated for this study are associated with NCBI Gene Expression Omnibus
accession **GSE223032**. Do not place secure reviewer tokens in this public code
repository.

The human oligodendrocyte analyses reuse published single-nucleus RNA-sequencing
datasets. Add the exact accession identifier, source publication, downloaded
file/version, and access date for each reused dataset before public release.

Large Seurat objects, raw count matrices, SCENIC databases, and derived result
directories are intentionally excluded from Git. Their repository locations and
identifiers should be recorded here once the final deposits are public.

Expected local inputs include:

- the 108,166-cell mouse singlet-only Seurat object after sample-wise
  scDblFinder filtering;
- the integrated human single-nucleus Seurat object;
- the human oligodendrocyte resolution-0.4 object;
- human and mouse cisTarget ranking databases listed in the SCENIC script; and
- machine-readable source tables for the manuscript panels listed in
  `docs/figure_script_map.tsv`.
