# Reproducible workflow

Run every command from the repository root after creating `config/paths.R`.
The scripts retain the analysis logic of the archived working versions while
using a central path configuration.

## 0. Raw data and quality control

Status: **incomplete in the local archive**. Recover the raw import, threshold
QC, and sample-wise scDblFinder scripts from the server. Do not claim an
end-to-end raw-data workflow until these scripts and their exact outputs are
present.

## 1. Mouse major-cell-type annotation

```bash
Rscript scripts/01_mouse_annotation/01_annotate_major_cell_types.R
```

Input: clustered mouse Seurat object. Output: annotated mouse Seurat object,
cell-type tables, and annotation plots.

## 2. Mouse oligodendrocyte analysis

```bash
Rscript scripts/02_mouse_oligodendrocytes/01_subcluster_and_annotate_mouse_ol.R
Rscript scripts/02_mouse_oligodendrocytes/02_subtype_markers_and_enrichment.R
```

These scripts extract oligodendrocytes, run sample-aware Harmony integration,
assign Oligo1/Oligo2/Oligo3/A4-OL labels, and calculate subtype markers and
enrichment results.

## 3. Human major-cell-type and oligodendrocyte analyses

```bash
Rscript scripts/03_human_annotation/01_annotate_major_cell_types.R
Rscript scripts/04_human_oligodendrocytes/01_resolution_sweep.R
Rscript scripts/04_human_oligodendrocytes/02_assign_a4_ol_subtypes.R
Rscript scripts/04_human_oligodendrocytes/03_a4_ol_markers_and_enrichment.R
```

Human cell identities and oligodendrocyte states are assigned from joint marker
expression patterns rather than a single marker gene.

## 4. Cross-species UCell validation

```bash
Rscript scripts/05_cross_species_ucell/01_mouse_signature_to_human_ucell.R
Rscript scripts/05_cross_species_ucell/02_plot_human_ucell_comparisons.R
```

The first script derives the top 200 positive mouse A4-OL markers using a
Wilcoxon rank-sum test, converts mouse markers to human orthologs with
`babelgene`, calculates UCell scores in human oligodendrocytes, and writes an
ortholog-mapping audit table.

## 5. SCENIC

```bash
Rscript scripts/06_scenic/01_run_scenic_human_mouse.R
Rscript scripts/06_scenic/02_plot_egr1_regulon_volcano.R
Rscript scripts/06_scenic/03_egr1_downstream_summary.R
```

SCENIC requires species-specific cisTarget ranking databases. Configure their
parent directory using `A4OL_CISTARGET_ROOT`.

## 6. Mouse validation and robustness

```bash
Rscript scripts/07_mouse_validation/01_mouse_a4_ucell_by_genotype.R
Rscript scripts/07_mouse_validation/02_leave_one_mouse_out_robustness.R
```

The leave-one-mouse-out workflow reruns normalization, PCA, sample-identity
Harmony integration, clustering, signature derivation, and UCell scoring after
excluding each mouse in turn. APOE genotype is not supplied as a Harmony
correction variable.

## 7. Final figure rendering

Run only the plotting scripts required by the final manuscript panels. Their
inputs are outputs of the SCENIC and UCell stages. The current main and Extended
Data panel identifiers are recorded in `docs/figure_script_map.tsv`.
