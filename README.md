# Cross-species analysis of A4-like oligodendrocytes and APOE genotype

This repository contains the curated analysis code used to identify and
characterize A4-like oligodendrocytes (A4-OLs) in mouse and human single-cell or
single-nucleus RNA-sequencing data. It includes cell-type annotation,
oligodendrocyte subclustering, differential-expression and enrichment analyses,
mouse-to-human ortholog mapping, UCell scoring, SCENIC regulon analysis, and
leave-one-mouse-out robustness analyses.

The public repository is organized by analysis stage. Superseded, exploratory,
debugging, and plotting-repair scripts from the working directory are not part
of the reproducible main workflow. The original-to-public filename mapping is
recorded in `docs/original_script_map.tsv`.

## Repository structure

```text
.
├── R/                          Shared configuration loader
├── config/                     Portable path configuration template
├── metadata/                   Sample identifiers and QC counts
├── environment/                Dependency and session-information tools
├── docs/                       Workflow, file maps, and release checklist
├── scripts/
│   ├── 00_raw_data_and_qc/     Status of upstream raw-data/QC scripts
│   ├── 01_mouse_annotation/
│   ├── 02_mouse_oligodendrocytes/
│   ├── 03_human_annotation/
│   ├── 04_human_oligodendrocytes/
│   ├── 05_cross_species_ucell/
│   ├── 06_scenic/
│   ├── 07_mouse_validation/
│   └── 08_figure_generation/
└── source_data/                Instructions for figure-level source data
```

## Configuration

All commands below assume that the current working directory is the repository
root.

1. Copy `config/paths.example.R` to `config/paths.R`.
2. Edit the private server roots in `config/paths.R`, or provide the equivalent
   environment variables.
3. Do not commit `config/paths.R`; it is intentionally ignored by Git.

Example:

```bash
cp config/paths.example.R config/paths.R
Rscript scripts/05_cross_species_ucell/01_mouse_signature_to_human_ucell.R
```

The main configuration variables are:

- `A4OL_SERVER_ROOT`: top-level server user/project root.
- `A4OL_PROJECT_ROOT`: root of the A4-OL sequencing project.
- `A4OL_RESULTS_ROOT`: root containing processed analysis outputs.
- `A4OL_FIGURE_ROOT`: root containing figure-analysis outputs.
- `A4OL_CISTARGET_ROOT`: directory containing human and mouse cisTarget files.
- `A4OL_FONT_FILE`: path to the TrueType font used to recreate figure typography.

## Analysis order

The intended order is:

1. Raw-data processing, cell QC, and sample-wise doublet detection.
2. Mouse and human major-cell-type annotation.
3. Mouse and human oligodendrocyte subclustering and annotation.
4. Mouse A4-OL marker derivation and ortholog-based mapping to human genes.
5. Human UCell validation.
6. Human and mouse SCENIC analyses.
7. Mouse genotype and leave-one-mouse-out robustness analyses.
8. Final figure rendering.

See `docs/workflow.md` for script-level dependencies and
`docs/figure_script_map.tsv` for output-to-script mapping.

## Manuscript panels covered

The panel map follows the current manuscript:

- Main Fig. 2a,c and Fig. 3b,c;
- Extended Data Fig. 5a-c and Fig. 6a,b;
- Extended Data Fig. 7a-h;
- Extended Data Fig. 8a,d-f.

Every listed panel has a corresponding analysis and/or final rendering script
in `docs/figure_script_map.tsv`.

## Mouse-to-human gene mapping

The cross-species UCell workflow uses `babelgene::orthologs()` to map mouse gene
symbols to human orthologs. It does not infer human symbols by converting mouse
symbols to uppercase. The script writes a mapping audit table containing mapped,
unmapped, one-to-many, and expression-matched genes.

## Quality-control status

The locally archived code begins after upstream raw-data processing. The
original raw-data import, threshold-based cell QC, and sample-wise scDblFinder
scripts must still be recovered from the analysis server before public release.
The verified doublet-filtering counts are provided in
`metadata/mouse_sample_key_and_doublet_qc.tsv` and summarized in
`scripts/00_raw_data_and_qc/README.md`.

## Software environment

Package names used by the repository are listed in
`environment/package_requirements.tsv`. Run
`Rscript environment/capture_session_info.R` in the final server environment to
record exact package versions and session information. A release-specific
`renv.lock` should then be generated and committed; no package versions have
been guessed in this repository.

## Data and code availability

The generated mouse single-nucleus RNA-sequencing and spatial transcriptomics
data are associated with GEO accession `GSE223032`. Processed objects and figure
source tables are not committed to Git. Reused human dataset identifiers and the
final code archive DOI still need to be completed. See `docs/data_access.md` and
`docs/release_checklist.md`; do not place secure reviewer tokens or temporary
cloud links in the public repository.

## License and citation

A software license and final `CITATION.cff` require author/institutional approval
and have therefore not been invented here. The required decisions are listed in
`docs/release_checklist.md`.

Git initialization, first push, version tagging, and archival steps are listed
in `docs/github_release.md`.
