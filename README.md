# Cross-species analysis of A4-like oligodendrocytes and APOE genotype

This repository contains the curated R code used to identify and characterize
A4-like oligodendrocytes (A4-OLs) in mouse and human single-cell or
single-nucleus RNA-sequencing data. The workflow covers raw mouse 10x input,
quality control, sample-wise doublet detection, major-cell-type annotation,
oligodendrocyte subclustering, A4-OL marker and signature derivation,
mouse-to-human ortholog mapping, human UCell validation, SCENIC regulon
analysis, APOE-group comparisons, leave-one-mouse-out (LOMO) robustness
analysis, enrichment, and final figure rendering.

The public repository is organized by analysis stage. Superseded,
exploratory, debugging, repair-only, and unrelated scripts from the working
directories are not part of the reproducible manuscript workflow. The
original-to-public filename mapping and the supplemental-code audit are in
`docs/original_script_map.tsv` and `docs/supplemental_code_audit.md`.

## Repository structure

```text
.
├── R/                          Shared configuration and ortholog helpers
├── config/                     Portable path configuration template
├── demo/                       Small deterministic cross-species demo
├── metadata/                   Mouse labels and verified doublet counts
├── environment/                Dependencies, tested platform, installation
├── docs/                       Workflow, figure map, audit, checklist
├── scripts/
│   ├── 00_raw_data_and_qc/     10x import, threshold QC, scDblFinder
│   ├── 01_mouse_annotation/
│   ├── 02_mouse_oligodendrocytes/
│   ├── 03_human_annotation/
│   ├── 04_human_oligodendrocytes/
│   ├── 05_cross_species_ucell/
│   ├── 06_scenic/
│   ├── 07_mouse_validation/
│   └── 08_figure_generation/
└── source_data/                Figure-level source-data guidance
```

## Tested system and dependencies

The available CFFF server session record documents:

- Ubuntu 22.04.3 LTS on x86_64 Linux;
- R 4.5.2;
- Seurat 5.4.0 and SeuratObject 5.3.0;
- 8 CPU cores and 250 GB RAM for full mouse preprocessing.

The scripts use portable `file.path()` configuration. Windows 10/11 is also
supported: scDblFinder falls back to a serial backend, SCENIC uses
`doParallel`, and no Unix `find` command is required. Full Seurat and SCENIC
analyses were run on the CFFF Linux server; a complete end-to-end Windows rerun
has not been claimed.

All direct dependencies are listed in
`environment/package_requirements.tsv`, with the complete exact versions used
in the tested environment dated 28 February 2026. CRAN packages are fixed to
that day's Posit CRAN snapshot, Bioconductor packages are fixed to
Bioconductor 3.22 revisions at or before that date, and SCENIC is fixed to its
corresponding official GitHub commit. The additional
`environment/cfff_linux_session_info.txt` file records a later CFFF server
session and is not the version authority for the dated tested environment.

## Installation

Install R 4.5.2 or recreate the documented CFFF environment, then run. The
installer selects the 2026-02-28 CRAN snapshot and Bioconductor 3.22 and checks
the resulting direct-package versions against the declared table:

```bash
Rscript environment/install_dependencies.R
```

Typical dependency installation is estimated at 15-45 minutes on a broadband
connection, excluding R and the large species-specific cisTarget database
files. Installation time varies by operating system and whether binary
packages are available. The code does not require a GPU.

## Configuration

All commands assume that the current working directory is the repository root.

1. Copy `config/paths.example.R` to `config/paths.R`.
2. Edit the private data and output roots, or set the equivalent environment
   variables.
3. Do not commit `config/paths.R`; it is intentionally ignored by Git.

Linux/macOS:

```bash
cp config/paths.example.R config/paths.R
```

Windows PowerShell:

```powershell
Copy-Item config/paths.example.R config/paths.R
```

The main configuration variables are:

- `A4OL_SERVER_ROOT`: top-level project/server root;
- `A4OL_PROJECT_ROOT`: A4-OL sequencing project root;
- `A4OL_RESULTS_ROOT`: processed analysis output root;
- `A4OL_FIGURE_ROOT`: figure-analysis output root;
- `A4OL_CISTARGET_ROOT`: human and mouse cisTarget parent directory;
- `A4OL_FONT_FILE`: TrueType font used for submitted figures;
- `A4OL_MOUSE_RAW_ROOT`: Cell Ranger output parent directory;
- `A4OL_MOUSE_QC_ROOT`: threshold-QC and scDblFinder output directory;
- `A4OL_N_CORES`: worker count for supported parallel steps.

## Quick demo

The small demo uses base R and does not require the full expression objects:

```bash
Rscript demo/run_demo.R
```

It should finish in under one minute on a normal desktop and create
`demo/output/ortholog_mapping_audit.tsv` and
`demo/output/human_signature.txt`. Expected files and instructions are in
`demo/expected_output/` and `demo/README.md`. The demo uses a frozen
ortholog fixture for determinism; the full analysis uses `babelgene` and
UCell.

## Full analysis order

After configuration, run the modules in this order:

1. raw mouse 10x import, threshold QC, and sample-wise scDblFinder;
2. mouse and human major-cell-type annotation;
3. mouse and human oligodendrocyte subclustering and annotation;
4. mouse A4-OL marker derivation and ortholog-based mapping to human genes;
5. human UCell validation;
6. human and mouse SCENIC analyses;
7. mouse APOE-group and LOMO robustness analyses; and
8. final figure rendering.

Commands, inputs, and outputs are listed in `docs/workflow.md`.
`docs/figure_script_map.tsv` maps every reported panel to its analysis and
rendering scripts.

## Manuscript panels covered

- Main Fig. 2a,c and Fig. 3b,c;
- Extended Data Fig. 5a-c and Fig. 6a,b;
- Extended Data Fig. 7a-h;
- Extended Data Fig. 8a,d-f.

## Mouse quality control and doublet removal

The raw-data script applies 700 < nFeature_RNA < 5,000,
500 < nCount_RNA < 30,000, mitochondrial fraction < 5%, and
log10(features)/log10(counts) > 0.75. Harmony is grouped by `sample`;
APOE genotype is not a correction variable.

The archived post-threshold object contained 120,057 cells. Sample-wise
scDblFinder removed 11,891 doublets (9.90%) and retained 108,166 singlets for
downstream analysis. Per-mouse counts are in
`metadata/mouse_sample_key_and_doublet_qc.tsv`.

## Mouse-to-human mapping

The full cross-species workflow maps mouse symbols to human orthologs with
`babelgene::orthologs()`. It does not infer orthologs by converting mouse
symbols to uppercase. The workflow writes an audit table containing mapped,
unmapped, one-to-many, and expression-matched genes before UCell scoring.

## Data, runtime, and hardware

Large Seurat objects, raw count matrices, cisTarget databases, and derived
results are intentionally excluded from Git. Configure their local locations
in `config/paths.R`. The mouse data are associated with GEO accession
`GSE223032`; human dataset identifiers and final source-data deposits remain
listed as release actions in `docs/data_access.md`.

The demo needs only a normal desktop. Full mouse preprocessing and SCENIC are
HPC-scale analyses; the archived preprocessing run used 8 CPU cores and
250 GB RAM. Full runtimes were not recovered and are therefore not invented.

## Release status

The source code, installation instructions, tested Linux/R versions, hardware
description, small demo, expected demo output, demo runtime, and instructions
for applying the code to configured data are present. A complete package
lockfile, end-to-end clean rerun, public figure source data, final human dataset
identifiers, citation metadata, versioned release, and archive DOI are still
outstanding. See `docs/software_submission_checklist.md` and
`docs/release_checklist.md`.

## License and citation

The code is released under the MIT License; see `LICENSE`. A final
`CITATION.cff` will be added after the manuscript title, author list, ORCIDs,
version, and archival DOI are finalized.
