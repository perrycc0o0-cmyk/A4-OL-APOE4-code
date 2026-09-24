# Shared configuration loader for all analysis scripts.
# Run scripts from the repository root, or set A4OL_CONFIG_FILE to an
# absolute path pointing to a private configuration file.

config_file <- Sys.getenv(
  "A4OL_CONFIG_FILE",
  unset = file.path("config", "paths.R")
)

if (!file.exists(config_file)) {
  stop(
    "Configuration file not found: ", config_file, "\n",
    "Copy config/paths.example.R to config/paths.R and edit it, or set ",
    "A4OL_CONFIG_FILE to an absolute path."
  )
}

sys.source(config_file, envir = .GlobalEnv)

required_config_objects <- c(
  "A4OL_SERVER_ROOT",
  "A4OL_PROJECT_ROOT",
  "A4OL_RESULTS_ROOT",
  "A4OL_FIGURE_ROOT",
  "A4OL_CISTARGET_ROOT",
  "A4OL_FONT_FILE"
)

missing_config_objects <- required_config_objects[
  !vapply(required_config_objects, exists, logical(1), inherits = TRUE)
]

if (length(missing_config_objects) > 0) {
  stop(
    "Missing configuration object(s): ",
    paste(missing_config_objects, collapse = ", ")
  )
}

if (!nzchar(A4OL_SERVER_ROOT)) {
  stop("A4OL_SERVER_ROOT is empty. Edit config/paths.R before running analyses.")
}

