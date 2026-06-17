options(stringsAsFactors = FALSE)

install_and_load <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

install_and_load("ggplot2")
install_and_load("ggVennDiagram")

args_all <- commandArgs(trailingOnly = FALSE)
file_arg_idx <- grep("^--file=", args_all)
if (length(file_arg_idx) > 0) {
  script_path <- sub("^--file=", "", args_all[file_arg_idx[1]])
  base_dir <- dirname(normalizePath(script_path))
} else {
  base_dir <- getwd()
}

input_dir <- file.path(dirname(base_dir), "converted_gene_id")
out_dir <- base_dir

input_files <- list(
  GSE_40 = "mapped_GSE_40_02_resLFC_Treatment_colon_cancer_vs_normal_biomart.csv",
  GSE_45 = "mapped_GSE_45_OR_02_resLFC_condition_tumor_vs_normal_biomart.csv",
  GSE_60 = "mapped_GSE_60_resLFC_raw_wol_Treatment_tumor_vs_normal_biomart.csv",
  GSE_79 = "mapped_GSE_79_02_resLFC_Treatment_tumor_vs_normal_biomart.csv",
  TCGA_COAD = "mapped_tcga_resLFC_Treatment_Tumor_vs_Normal_biomart.csv",
  TCGA_READ = "mapped_tcga_read_resLFC_Treatment_Tumor_vs_Normal_biomart.csv"
)

load_gene_sets <- function() {
  gene_sets <- list()

  for (nm in names(input_files)) {
    path <- file.path(input_dir, input_files[[nm]])
    if (!file.exists(path)) {
      stop(paste("Missing input file:", path))
    }

    df <- read.csv(path, check.names = FALSE)
    if (!("gene_id" %in% colnames(df))) {
      stop(paste("Column 'gene_id' not found in", path))
    }

    genes <- as.character(df$gene_id)
    genes <- trimws(genes)
    genes <- genes[!is.na(genes)]
    genes <- genes[genes != ""]
    genes <- genes[toupper(genes) != "NA"]
    gene_sets[[nm]] <- unique(genes)
  }

  gene_sets
}

exact_intersection_rows <- function(gene_sets) {
  names_vec <- names(gene_sets)
  rows <- list()
  idx <- 1

  for (k in seq_along(names_vec)) {
    combs <- combn(names_vec, k, simplify = FALSE)

    for (subset_names in combs) {
      included <- subset_names
      excluded <- setdiff(names_vec, included)

      overlap <- Reduce(intersect, gene_sets[included])
      if (length(excluded) > 0) {
        overlap <- setdiff(overlap, Reduce(union, gene_sets[excluded]))
      }

      label <- paste(included, collapse = " & ")
      if (length(overlap) == 0) {
        rows[[idx]] <- data.frame(
          intersection = label,
          included_sets = paste(included, collapse = "|"),
          n_sets = length(included),
          gene_id = "",
          stringsAsFactors = FALSE
        )
        idx <- idx + 1
      } else {
        overlap <- sort(overlap)
        rows[[idx]] <- data.frame(
          intersection = rep(label, length(overlap)),
          included_sets = rep(paste(included, collapse = "|"), length(overlap)),
          n_sets = rep(length(included), length(overlap)),
          gene_id = overlap,
          stringsAsFactors = FALSE
        )
        idx <- idx + 1
      }
    }
  }

  do.call(rbind, rows)
}

exact_intersection_counts <- function(gene_sets) {
  names_vec <- names(gene_sets)
  rows <- list()
  idx <- 1

  for (k in seq_along(names_vec)) {
    combs <- combn(names_vec, k, simplify = FALSE)

    for (subset_names in combs) {
      included <- subset_names
      excluded <- setdiff(names_vec, included)
      overlap <- Reduce(intersect, gene_sets[included])
      if (length(excluded) > 0) {
        overlap <- setdiff(overlap, Reduce(union, gene_sets[excluded]))
      }

      rows[[idx]] <- data.frame(
        intersection = paste(included, collapse = " & "),
        included_sets = paste(included, collapse = "|"),
        n_sets = length(included),
        count = length(overlap),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
    }
  }

  out <- do.call(rbind, rows)
  out <- out[order(out$n_sets, out$intersection), ]
  rownames(out) <- NULL
  out
}

gene_sets <- load_gene_sets()
cat("Loaded gene_id sets (R):\n")
for (nm in names(gene_sets)) {
  cat(sprintf("  %s: %d\n", nm, length(gene_sets[[nm]])))
}

venn_plot <- ggVennDiagram(
  gene_sets,
  label = "count",
  label_geom = "text",
  label_alpha = 1,
  label_color = "black",
  label_size = 4,
  edge_size = 0.8,
  set_name_fill = NA,
  set_name_color = "black",
  set_name_size = 4
) +
  ggplot2::scale_fill_gradient(low = "#e7f4f8", high = "#0072b2") +
  ggplot2::ggtitle("Gene Overlap Across 6 Datasets (mapped gene_id)")

png_path <- file.path(out_dir, "venn_6sets_gene_id_R.png")
pdf_path <- file.path(out_dir, "venn_6sets_gene_id_R.pdf")
ggplot2::ggsave(png_path, plot = venn_plot, width = 10, height = 8)
ggplot2::ggsave(pdf_path, plot = venn_plot, width = 10, height = 8)

cat("Saved R Venn PNG:", png_path, "\n")
cat("Saved R Venn PDF:", pdf_path, "\n")

intersections_df <- exact_intersection_rows(gene_sets)
intersections_path <- file.path(out_dir, "all_possible_exact_intersections_gene_lists_R.csv")
write.csv(intersections_df, intersections_path, row.names = FALSE)
cat("Saved intersections CSV (R):", intersections_path, "\n")

counts_df <- exact_intersection_counts(gene_sets)
counts_path <- file.path(out_dir, "exact_intersection_counts_R.csv")
write.csv(counts_df, counts_path, row.names = FALSE)
cat("Saved counts CSV (R):", counts_path, "\n")
