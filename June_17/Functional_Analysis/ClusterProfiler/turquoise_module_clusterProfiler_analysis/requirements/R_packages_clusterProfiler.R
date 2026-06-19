file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[[1]]) else file.path(getwd(), "requirements", "R_packages_clusterProfiler.R")
analysis_dir <- dirname(dirname(normalizePath(script_path, mustWork = FALSE)))
local_lib <- file.path(analysis_dir, "r_lib")
dir.create(local_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(local_lib, .libPaths()))

source(file.path(analysis_dir, "scripts", "install_clusterProfiler_packages.R"))
