args <- commandArgs(trailingOnly = TRUE)
analysis_dir <- if (length(args) >= 1) normalizePath(args[[1]], mustWork = FALSE) else normalizePath(file.path(getwd(), ".."), mustWork = FALSE)
local_lib <- file.path(analysis_dir, "r_lib")
dir.create(local_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(local_lib, .libPaths()))

cran_packages <- c(
  "readr",
  "dplyr",
  "tidyr",
  "stringr",
  "ggplot2",
  "forcats",
  "tibble",
  "purrr",
  "ggupset"
)

missing_cran <- setdiff(cran_packages, rownames(installed.packages()))
if (length(missing_cran) > 0) {
  install.packages(missing_cran, lib = local_lib, repos = "https://cloud.r-project.org")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = local_lib, repos = "https://cloud.r-project.org")
}

bioc_packages <- c(
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot",
  "ReactomePA",
  "DOSE",
  "AnnotationDbi"
)

missing_bioc <- setdiff(bioc_packages, rownames(installed.packages()))
if (length(missing_bioc) > 0) {
  BiocManager::install(missing_bioc, lib = local_lib, ask = FALSE, update = FALSE)
}

installed <- data.frame(
  package = c(cran_packages, bioc_packages),
  installed = vapply(c(cran_packages, bioc_packages), requireNamespace, logical(1), quietly = TRUE)
)
write.csv(installed, file.path(analysis_dir, "requirements", "installed_package_status.csv"), row.names = FALSE)
print(installed)
