#!/usr/bin/env Rscript

dir.create("environment", recursive = TRUE, showWarnings = FALSE)

sink(file.path("environment", "sessionInfo.txt"))
print(sessionInfo())
sink()

capture.output(
  list(
    R = R.version.string,
    platform = R.version$platform,
    os = Sys.info(),
    timezone = Sys.timezone()
  ),
  file = file.path("environment", "platform_info.txt")
)

installed <- as.data.frame(installed.packages(), stringsAsFactors = FALSE)
installed <- installed[, c("Package", "Version", "LibPath", "Priority")]
write.csv(
  installed,
  file.path("environment", "installed_packages.csv"),
  row.names = FALSE
)

requirements <- read.delim(
  file.path("environment", "package_requirements.tsv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
declared <- merge(
  requirements,
  installed[, c("Package", "Version")],
  by.x = "package",
  by.y = "Package",
  all.x = TRUE,
  sort = FALSE
)
write.table(
  declared,
  file.path("environment", "declared_package_versions.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message(
  "Wrote sessionInfo.txt, platform_info.txt, installed_packages.csv, and ",
  "declared_package_versions.tsv"
)
