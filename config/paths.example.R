# Copy this file to config/paths.R and edit the values for your server.
# config/paths.R is ignored by Git so private server locations are not shared.

A4OL_SERVER_ROOT <- Sys.getenv(
  "A4OL_SERVER_ROOT",
  unset = "/path/to/server/user/root"
)

A4OL_PROJECT_ROOT <- Sys.getenv(
  "A4OL_PROJECT_ROOT",
  unset = file.path(A4OL_SERVER_ROOT, "scRNAseq", "qin-snRNAseq")
)

A4OL_RESULTS_ROOT <- Sys.getenv(
  "A4OL_RESULTS_ROOT",
  unset = file.path(A4OL_PROJECT_ROOT, "result-qin")
)

A4OL_FIGURE_ROOT <- Sys.getenv(
  "A4OL_FIGURE_ROOT",
  unset = file.path(A4OL_RESULTS_ROOT, "Figure")
)

A4OL_CISTARGET_ROOT <- Sys.getenv(
  "A4OL_CISTARGET_ROOT",
  unset = file.path(A4OL_SERVER_ROOT, "database", "cisTarget")
)

# Set this to the font file used for the submitted figures.
A4OL_FONT_FILE <- Sys.getenv("A4OL_FONT_FILE", unset = "")

# Directory containing one Cell Ranger output folder per mouse, named
# <internal_sample_id>_outs/filtered_feature_bc_matrix.h5.
A4OL_MOUSE_RAW_ROOT <- Sys.getenv(
  "A4OL_MOUSE_RAW_ROOT",
  unset = file.path(A4OL_PROJECT_ROOT, "data", "20251211")
)

# Output directory for raw import, threshold QC, and scDblFinder reports.
A4OL_MOUSE_QC_ROOT <- Sys.getenv(
  "A4OL_MOUSE_QC_ROOT",
  unset = file.path(A4OL_RESULTS_ROOT, "raw_data_and_qc")
)

# Number of workers used by scripts that support parallel execution.
A4OL_N_CORES <- suppressWarnings(as.integer(Sys.getenv("A4OL_N_CORES", "8")))
if (is.na(A4OL_N_CORES) || A4OL_N_CORES < 1L) A4OL_N_CORES <- 1L
