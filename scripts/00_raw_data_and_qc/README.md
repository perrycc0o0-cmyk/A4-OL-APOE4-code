# Raw-data processing and quality control

The upstream workflow was reconstructed from the archived CFFF server scripts
in `analysis-refine` and `Afigureoutput+`, then made portable by moving private
paths into `config/paths.R`. Run from the repository root:

```bash
Rscript scripts/00_raw_data_and_qc/01_build_mouse_seurat_from_10x.R
Rscript scripts/00_raw_data_and_qc/02_samplewise_scdblfinder.R
```

The first script reads Cell Ranger `filtered_feature_bc_matrix.h5` files,
constructs one Seurat object per mouse, applies the archived thresholds, merges
the retained cells, and runs PCA, sample-identity Harmony, clustering, and UMAP.
APOE genotype is metadata only and is not a Harmony correction variable.

The applied cell filters use strict inequalities:

- 700 < detected features < 5,000;
- 500 < RNA counts < 30,000;
- mitochondrial fraction < 5%; and
- log10(features)/log10(counts) > 0.75.

`CreateSeuratObject` additionally uses `min.cells = 3` and
`min.features = 300`. The script writes the thresholds and per-sample counts
before and after threshold QC to machine-readable TSV files.

The verified sample-wise scDblFinder summary is:

- input object: 120,057 cells;
- retained singlets: 108,166 cells;
- removed doublets: 11,891 cells (9.90%);
- downstream analyses used the 108,166-cell singlet-only object.

Per-sample counts are recorded in
`metadata/mouse_sample_key_and_doublet_qc.tsv`.

`02_samplewise_scdblfinder.R` runs `scDblFinder` with `samples =
"doublet_qc_sample"`, `clusters = TRUE`, `dbr = NULL`, and
`multiSampleMode = "split"`. Linux uses `MulticoreParam`; Windows uses the
portable serial backend. The singlet-only output is also saved as
`seurat_clustered.rds`, which is the input to major-cell-type annotation.
