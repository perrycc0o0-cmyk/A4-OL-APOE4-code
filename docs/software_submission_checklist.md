# Nature Research code and software checklist status

This audit uses the repository contents as evidence. A box is checked only when
the requested item is actually present.

## Code and README items

- [x] Source code is provided.
- [x] A small simulated/frozen demo dataset is provided.
- [x] All direct dependencies and operating systems are listed with complete
      version numbers. Ubuntu 22.04.3 LTS and R 4.5.2 are documented;
      `environment/package_requirements.tsv` provides exact package versions
      for the tested environment dated 2026-02-28.
- [x] Tested software/platform versions are stated: CFFF Ubuntu 22.04.3 LTS,
      R 4.5.2, Seurat 5.4.0, and SeuratObject 5.3.0. Windows 10/11 is described
      as compatible, not as a completed end-to-end validation.
- [x] Non-standard hardware requirements are stated. Full preprocessing used
      8 CPU cores and 250 GB RAM; no GPU is required.
- [x] Installation instructions are provided.
- [x] Typical dependency-installation time is provided as an estimate.
- [x] Instructions to run the software on a small demo are provided.
- [x] Expected demo output is provided.
- [x] Expected demo runtime is provided.
- [x] Instructions for applying the code to configured user data are provided.
- [ ] Complete instructions and deposited inputs sufficient to reproduce every
      manuscript figure end to end are not yet available. The remaining items
      are the source-data deposit, human dataset identifiers, and clean
      full-workflow rerun.

## Additional form fields

- Software licence: **MIT License**. The approved licence text is provided in
  the repository-level `LICENSE` file.
- Code repository:
  https://github.com/perrycc0o0-cmyk/A4-OL-APOE4-code
- Functionality described: select **Elsewhere** and enter
  “GitHub README, docs/workflow.md, and docs/figure_script_map.tsv”.

## Functionality to list in the form

“Raw 10x snRNA-seq import and quality control; sample-wise scDblFinder doublet
removal; mouse and human cell-type annotation; oligodendrocyte subclustering;
A4-OL marker and signature derivation; mouse-to-human ortholog mapping;
cross-species UCell validation; SCENIC regulon analysis; EGR1/Egr1 downstream
analysis; APOE genotype comparisons; leave-one-mouse-out robustness analysis;
GO/KEGG enrichment; and manuscript figure generation.”
