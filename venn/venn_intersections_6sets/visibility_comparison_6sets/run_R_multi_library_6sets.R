options(stringsAsFactors = FALSE)
options(repos = c(CRAN = "https://cloud.r-project.org"))

install_and_load <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

install_and_load("ggplot2")
install_and_load("ggVennDiagram")
install_and_load("eulerr")
install_and_load("UpSetR")

args_all <- commandArgs(trailingOnly = FALSE)
file_arg_idx <- grep("^--file=", args_all)
if (length(file_arg_idx) > 0) {
  script_path <- sub("^--file=", "", args_all[file_arg_idx[1]])
  base_dir <- dirname(normalizePath(script_path))
} else {
  base_dir <- getwd()
}

input_dir <- file.path(dirname(base_dir), "..", "converted_gene_id")
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
    p <- file.path(input_dir, input_files[[nm]])
    df <- read.csv(p, check.names = FALSE)
    genes <- as.character(df$gene_id)
    genes <- trimws(genes)
    genes <- genes[!is.na(genes)]
    genes <- genes[genes != ""]
    genes <- genes[toupper(genes) != "NA"]
    gene_sets[[nm]] <- unique(genes)
  }
  gene_sets
}

gene_sets <- load_gene_sets()

# 1) ggVennDiagram
p1 <- ggVennDiagram(
  gene_sets,
  label = "count",
  label_geom = "label",
  label_alpha = 1,
  label_fill = "black",
  label_color = "white",
  label_size = 3.0,
  edge_size = 0.8
) +
  ggplot2::ggtitle("R: ggVennDiagram (6 sets)")

ggsave(file.path(out_dir, "R_ggVennDiagram_6sets.png"), plot = p1, width = 10, height = 8)
ggsave(file.path(out_dir, "R_ggVennDiagram_6sets.pdf"), plot = p1, width = 10, height = 8)

# 2) eulerr
fit <- eulerr::euler(gene_sets)

png(file.path(out_dir, "R_eulerr_6sets.png"), width = 1700, height = 1300, res = 150)
plot(
  fit,
  quantities = TRUE,
  labels = TRUE,
  fills = list(fill = c("#66c2a5", "#fc8d62", "#e78ac3", "#a6d854", "#e5c494", "#b3b3b3"), alpha = 0.45),
  main = "R: eulerr (6 sets)"
)
dev.off()

pdf(file.path(out_dir, "R_eulerr_6sets.pdf"), width = 12, height = 9)
plot(
  fit,
  quantities = TRUE,
  labels = TRUE,
  fills = list(fill = c("#66c2a5", "#fc8d62", "#e78ac3", "#a6d854", "#e5c494", "#b3b3b3"), alpha = 0.45),
  main = "R: eulerr (6 sets)"
)
dev.off()

# 3) UpSetR (intersection-focused visibility)
upset_input <- UpSetR::fromList(gene_sets)

png(file.path(out_dir, "R_UpSetR_6sets.png"), width = 1800, height = 1100, res = 160)
UpSetR::upset(
  upset_input,
  sets = names(gene_sets),
  keep.order = TRUE,
  order.by = "freq",
  nintersects = 50,
  mb.ratio = c(0.6, 0.4),
  point.size = 2.6,
  line.size = 0.8,
  text.scale = c(1.3, 1.2, 1.2, 1.0, 1.2, 1.2)
)
dev.off()

pdf(file.path(out_dir, "R_UpSetR_6sets.pdf"), width = 14, height = 8)
UpSetR::upset(
  upset_input,
  sets = names(gene_sets),
  keep.order = TRUE,
  order.by = "freq",
  nintersects = 50,
  mb.ratio = c(0.6, 0.4),
  point.size = 2.6,
  line.size = 0.8,
  text.scale = c(1.3, 1.2, 1.2, 1.0, 1.2, 1.2)
)
dev.off()

summary_path <- file.path(out_dir, "R_plot_summary.txt")
cat("R 6-set visibility comparison\n", file = summary_path)
cat("============================\n\n", file = summary_path, append = TRUE)
cat("Input set sizes:\n", file = summary_path, append = TRUE)
for (nm in names(gene_sets)) {
  cat(sprintf("- %s: %d\n", nm, length(gene_sets[[nm]])), file = summary_path, append = TRUE)
}
cat("\nGenerated plots:\n", file = summary_path, append = TRUE)
cat("- ggVennDiagram: R_ggVennDiagram_6sets.png/.pdf\n", file = summary_path, append = TRUE)
cat("- eulerr: R_eulerr_6sets.png/.pdf\n", file = summary_path, append = TRUE)
cat("- UpSetR: R_UpSetR_6sets.png/.pdf\n", file = summary_path, append = TRUE)

cat("R multi-library plots generated successfully.\n")
