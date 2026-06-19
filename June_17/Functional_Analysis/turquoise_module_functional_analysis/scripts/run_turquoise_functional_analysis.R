suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(gprofiler2)
  library(ggplot2)
  library(forcats)
})

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[[1]]) else "."
analysis_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
input_file <- file.path(analysis_dir, "inputs", "turquoise_module_genes.csv")
processed_dir <- file.path(analysis_dir, "results", "processed")
raw_dir <- file.path(analysis_dir, "results", "raw")
png_dir <- file.path(analysis_dir, "figures", "png")
pdf_dir <- file.path(analysis_dir, "figures", "pdf")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(png_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)

genes <- read_csv(input_file, show_col_types = FALSE) |>
  separate_rows(gene, sep = ";") |>
  mutate(gene = str_trim(gene)) |>
  filter(gene != "") |>
  distinct(gene) |>
  arrange(gene)

write_csv(genes, file.path(processed_dir, "turquoise_module_clean_genes_from_R.csv"))

gp <- gost(
  query = genes$gene,
  organism = "hsapiens",
  sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "WP"),
  correction_method = "g_SCS",
  user_threshold = 0.05,
  evcodes = TRUE
)

results <- gp$result |>
  arrange(p_value) |>
  mutate(
    negative_log10_p = -log10(p_value),
    intersection_genes = vapply(intersections, paste, collapse = ";", FUN.VALUE = character(1))
  )

write_csv(results, file.path(processed_dir, "gprofiler_enrichment_all_results_from_R.csv"))
write_csv(filter(results, p_value < 0.05), file.path(processed_dir, "gprofiler_enrichment_significant_from_R.csv"))

top_terms <- results |>
  filter(p_value < 0.05) |>
  slice_min(p_value, n = 25) |>
  mutate(label = paste0("[", source, "] ", term_name))

p <- ggplot(top_terms, aes(x = negative_log10_p, y = fct_reorder(label, negative_log10_p))) +
  geom_col(fill = "#4C78A8") +
  labs(
    title = "Turquoise module: top enriched functional terms",
    x = "-log10(g:SCS adjusted p-value)",
    y = NULL
  ) +
  theme_bw(base_size = 12)

ggsave(file.path(png_dir, "top25_enriched_terms_barplot_from_R.png"), p, width = 11, height = 8, dpi = 300)
ggsave(file.path(pdf_dir, "top25_enriched_terms_barplot_from_R.pdf"), p, width = 11, height = 8)
