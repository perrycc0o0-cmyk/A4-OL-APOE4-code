#!/usr/bin/env Rscript

dir.create("environment", recursive = TRUE, showWarnings = FALSE)

sink(file.path("environment", "sessionInfo.txt"))
print(sessionInfo())
sink()

installed <- as.data.frame(installed.packages(), stringsAsFactors = FALSE)
installed <- installed[, c("Package", "Version", "LibPath", "Priority")]
write.csv(
  installed,
  file.path("environment", "installed_packages.csv"),
  row.names = FALSE
)

message("Wrote environment/sessionInfo.txt and environment/installed_packages.csv")

