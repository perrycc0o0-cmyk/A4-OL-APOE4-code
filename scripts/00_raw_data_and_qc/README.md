# Raw-data processing and QC: scripts to recover before release

The available local archive does not contain the upstream scripts that import
raw single-cell/single-nucleus RNA-sequencing matrices, construct the initial
Seurat object, apply feature/count/mitochondrial thresholds, or run scDblFinder.
Those scripts must be copied from the analysis server into this directory before
the repository is declared complete.

The verified sample-wise scDblFinder summary is:

- input object: 120,057 cells;
- retained singlets: 108,166 cells;
- removed doublets: 11,891 cells (9.90%);
- downstream analyses used the 108,166-cell singlet-only object.

Per-sample counts are recorded in
`metadata/mouse_sample_key_and_doublet_qc.tsv`.

The recovered scripts should report, without exposing private filesystem paths:

1. raw input format and reference genome/annotation;
2. minimum and maximum detected-feature thresholds;
3. minimum and maximum count thresholds;
4. mitochondrial-read threshold;
5. whether thresholds differed by sample;
6. the exact scDblFinder call and package version;
7. cells/nuclei removed at each QC step; and
8. the filename of the 108,166-cell singlet-only output object.

