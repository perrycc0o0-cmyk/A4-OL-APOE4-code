# Small cross-species signature demo

This lightweight example exercises the table-level part of the cross-species
workflow without requiring the full Seurat objects or cisTarget databases. It
selects positive mouse A4-OL markers, joins a frozen mouse-human ortholog table,
checks which human orthologs occur in a small human feature list, and writes a
mapping audit plus the resulting human signature.

From the repository root, run:

```bash
Rscript demo/run_demo.R
```

Inputs are in `demo/input/`. Outputs are written to `demo/output/` and should
match the committed files in `demo/expected_output/`:

- `ortholog_mapping_audit.tsv`: six selected positive mouse markers, their
  frozen human orthologs, and presence/absence in the human feature list;
- `human_signature.txt`: `YWHAH`, `CAMK2N1`, `EGR1`, `FOS`, and `JUN`, in the
  marker-ranking order used by the demo.

Expected runtime is under one minute on a normal desktop after R is installed;
the script uses base R only. The frozen ortholog table makes the test
deterministic. The full analysis does not use this fixture: it calls
`babelgene::orthologs()` and then UCell on the complete human expression object.

To adapt the demo to another marker list, preserve the three required columns
in `mouse_a4_markers.tsv` (`mouse_symbol`, `avg_log2FC`, and `p_val_adj`) and
replace the ortholog and detected-feature files with the corresponding species
mapping and expression features.
