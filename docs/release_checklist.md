# Public-release checklist

## Blocking items

- [ ] Recover raw-data import, threshold QC, and sample-wise scDblFinder scripts
      from the server and place them in `scripts/00_raw_data_and_qc/`.
- [ ] Record the exact feature, count, and mitochondrial thresholds and the
      number of cells/nuclei removed at each QC step.
- [ ] Run every curated script from a clean output directory using the final
      server configuration.
- [ ] Run `Rscript environment/capture_session_info.R` on the server.
- [ ] Create and commit the final `renv.lock` from that verified environment.
- [ ] Recheck `docs/figure_script_map.tsv` if panel numbering changes during
      revision.
- [ ] Add figure-level machine-readable source data.
- [ ] Add raw/processed sequencing accession numbers and reviewer-access links.
- [ ] Confirm that no private human metadata or server credentials are tracked.
- [ ] Select a software license with the authors/institution (for example MIT or
      BSD-3-Clause) and add the approved `LICENSE` file.
- [ ] Add a valid `CITATION.cff` with the final title, authors, ORCIDs, version,
      repository URL, and release DOI.
- [ ] Create a versioned GitHub release and archive that release in Zenodo or an
      equivalent repository that issues a persistent DOI.

## Consistency checks

- [ ] README, Methods, Data Availability, and Code Availability use the same
      sample labels, software versions, data identifiers, and analysis order.
- [ ] Mouse-to-human mapping is ortholog based; no uppercase-only conversion
      remains in the public workflow.
- [ ] Harmony is grouped by sample identity; APOE genotype is not a correction
      variable.
- [ ] The repository contains only the final scripts needed to reproduce the
      reported results; exploratory scripts are retained outside the release.
