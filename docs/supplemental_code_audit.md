# Audit of supplemental working-code directories

Audit date: 24 September 2026.

## Scope

- `Afigureoutput+`: 103 files, including final, fixed, superseded, repair-only,
  test, log, archive, and unrelated analysis files.
- `analysis-refine`: 21 files covering early raw processing, integration,
  annotation, marker, and oligodendrocyte refinement work.

The working directories remain unchanged. Only curated copies needed for the
reported A4-OL manuscript workflow are included here.

## Files promoted into the public workflow

- `analysis-refine/process_harmony.R` supplied the raw 10x import, threshold
  QC, PCA, Harmony, clustering, and output logic for
  `scripts/00_raw_data_and_qc/01_build_mouse_seurat_from_10x.R`.
- `Afigureoutput+/QC_MOUSE_seurat_final_scDblFinder_BY_SAMPLE_FIXED.R`
  supplied the sample-wise doublet workflow for
  `scripts/00_raw_data_and_qc/02_samplewise_scdblfinder.R`.
- Final A4-OL signature, UCell, SCENIC, LOMO, annotation, enrichment, and panel
  scripts already curated in the repository remain mapped in
  `docs/original_script_map.tsv`.

## Corrections made during curation

- Removed private absolute CFFF and home-directory paths in favor of
  `config/paths.R`.
- Kept APOE genotype as metadata and used biological sample identity as the
  Harmony grouping variable.
- Preserved the archived threshold values and sample-wise scDblFinder settings.
- Replaced a Unix-only recursive `find` call with portable R file discovery.
- Added Windows-aware parallel backends: serial BiocParallel for scDblFinder
  and `doParallel` for SCENIC; Linux retains multicore backends.
- Added dependency installation, a captured CFFF platform record, a small demo,
  expected demo outputs, and checklist evidence.

## Files intentionally excluded

- Superseded variants with suffixes such as `FIXED`, `V2`, `V3`, or
  `V4` when a later curated public equivalent already exists.
- One-off redraw, formatting repair, package-install test, and debugging files
  that do not define independent analysis steps.
- Runtime logs, ZIP archives, generated outputs, and files containing private
  absolute paths.
- CHO/MPRA and TNR analyses, which are not mapped to the manuscript panels
  listed in `docs/figure_script_map.tsv`.
- The pseudobulk mouse A4-OL script, because that analysis is not part of the
  current Methods/panel scope.

## Remaining reproducibility gaps

- Exact tested versions of several Bioconductor/SCENIC/UCell dependencies must
  be captured from the final CFFF environment.
- The complete workflow has not yet been rerun from a clean output directory.
- Large input objects, human dataset identifiers, figure source-data tables,
  software licence, release tag, and archival DOI remain release actions.
