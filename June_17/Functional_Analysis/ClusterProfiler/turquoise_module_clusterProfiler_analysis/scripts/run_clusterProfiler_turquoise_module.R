args <- commandArgs(trailingOnly = TRUE)
analysis_dir <- if (length(args) >= 1) normalizePath(args[[1]], mustWork = FALSE) else normalizePath(file.path(getwd(), ".."), mustWork = FALSE)
local_lib <- file.path(analysis_dir, "r_lib")
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

required <- c(
  "readr", "dplyr", "tidyr", "stringr", "ggplot2", "forcats", "purrr", "ggupset",
  "clusterProfiler", "org.Hs.eg.db", "enrichplot", "ReactomePA", "DOSE", "AnnotationDbi"
)
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing, collapse = ", "),
    "\nRun scripts/install_clusterProfiler_packages.R first."
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(forcats)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ReactomePA)
  library(DOSE)
})

input_dir <- file.path(analysis_dir, "inputs")
raw_dir <- file.path(analysis_dir, "results", "raw")
processed_dir <- file.path(analysis_dir, "results", "processed")
png_dir <- file.path(analysis_dir, "figures", "png")
pdf_dir <- file.path(analysis_dir, "figures", "pdf")
report_dir <- file.path(analysis_dir, "report")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(png_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
unlink(file.path(raw_dir, "*_error.txt"))

input_file <- file.path(input_dir, "turquoise_module_genes.csv")
pvalue_cutoff <- 0.05
qvalue_cutoff <- 0.20
show_category <- 20

save_table <- function(object, filename) {
  table <- as.data.frame(object)
  readr::write_csv(table, file.path(processed_dir, filename))
  invisible(table)
}

save_plot <- function(plot, basename, width = 11, height = 8) {
  png_file <- file.path(png_dir, paste0(basename, ".png"))
  pdf_file <- file.path(pdf_dir, paste0(basename, ".pdf"))
  ggsave(png_file, plot = plot, width = width, height = height, dpi = 300, bg = "white")
  ggsave(pdf_file, plot = plot, width = width, height = height, bg = "white")
}

safe_plot <- function(expr, basename, width = 11, height = 8) {
  tryCatch({
    plot <- force(expr)
    save_plot(plot, basename, width, height)
    TRUE
  }, error = function(e) {
    writeLines(conditionMessage(e), file.path(raw_dir, paste0(basename, "_plot_error.txt")))
    FALSE
  })
}

clean_genes <- read_csv(input_file, show_col_types = FALSE) %>%
  tidyr::separate_rows(gene, sep = ";") %>%
  mutate(gene = str_trim(gene)) %>%
  filter(!is.na(gene), gene != "") %>%
  distinct(gene) %>%
  arrange(gene)

write_csv(clean_genes, file.path(processed_dir, "turquoise_module_clean_symbol_genes.csv"))
writeLines(clean_genes$gene, file.path(processed_dir, "turquoise_module_clean_symbol_genes.txt"))

symbol_to_entrez <- bitr(
  clean_genes$gene,
  fromType = "SYMBOL",
  toType = c("ENTREZID", "SYMBOL"),
  OrgDb = org.Hs.eg.db
) %>%
  distinct(SYMBOL, ENTREZID)

write_csv(symbol_to_entrez, file.path(processed_dir, "symbol_to_entrez_mapping.csv"))

unmapped <- clean_genes %>%
  anti_join(symbol_to_entrez, by = c("gene" = "SYMBOL"))
write_csv(unmapped, file.path(processed_dir, "unmapped_symbol_genes.csv"))

gene_entrez <- unique(symbol_to_entrez$ENTREZID)
writeLines(gene_entrez, file.path(processed_dir, "turquoise_module_entrez_genes.txt"))

run_metadata <- tibble::tibble(
  metric = c("input_symbol_genes", "mapped_entrez_ids", "unmapped_symbol_genes", "pvalue_cutoff", "qvalue_cutoff"),
  value = c(nrow(clean_genes), length(gene_entrez), nrow(unmapped), pvalue_cutoff, qvalue_cutoff)
)
write_csv(run_metadata, file.path(processed_dir, "run_metadata.csv"))

go_results <- list()
for (ontology in c("BP", "MF", "CC")) {
  ego <- enrichGO(
    gene = gene_entrez,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff = pvalue_cutoff,
    qvalueCutoff = qvalue_cutoff,
    readable = TRUE
  )
  go_results[[ontology]] <- ego
  table <- save_table(ego, paste0("enrichGO_", ontology, "_results.csv"))

  if (nrow(table) > 0) {
    safe_plot(dotplot(ego, showCategory = show_category, font.size = 10) + ggtitle(paste("GO", ontology, "enrichment - turquoise module")), paste0("dotplot_enrichGO_", ontology), 11, 8)
    safe_plot(barplot(ego, showCategory = show_category, font.size = 10) + ggtitle(paste("GO", ontology, "enrichment - turquoise module")), paste0("barplot_enrichGO_", ontology), 11, 8)
    safe_plot(upsetplot(ego, n = min(15, nrow(table))) + ggtitle(paste("GO", ontology, "gene overlap upset plot")), paste0("upsetplot_enrichGO_", ontology), 12, 8)
    ego_sim <- pairwise_termsim(ego)
    safe_plot(emapplot(ego_sim, showCategory = min(30, nrow(table))) + ggtitle(paste("GO", ontology, "enrichment map")), paste0("emapplot_enrichGO_", ontology), 12, 9)
  }
}

kegg_result <- tryCatch({
  enrichKEGG(
    gene = gene_entrez,
    organism = "hsa",
    keyType = "kegg",
    pAdjustMethod = "BH",
    pvalueCutoff = pvalue_cutoff,
    qvalueCutoff = qvalue_cutoff
  )
}, error = function(e) {
  writeLines(conditionMessage(e), file.path(raw_dir, "enrichKEGG_error.txt"))
  NULL
})

if (!is.null(kegg_result)) {
  kegg_readable <- setReadable(kegg_result, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  kegg_table <- save_table(kegg_readable, "enrichKEGG_results.csv")
  if (nrow(kegg_table) > 0) {
    safe_plot(dotplot(kegg_readable, showCategory = show_category, font.size = 10) + ggtitle("KEGG enrichment - turquoise module"), "dotplot_enrichKEGG", 11, 8)
    safe_plot(barplot(kegg_readable, showCategory = show_category, font.size = 10) + ggtitle("KEGG enrichment - turquoise module"), "barplot_enrichKEGG", 11, 8)
    safe_plot(upsetplot(kegg_readable, n = min(15, nrow(kegg_table))) + ggtitle("KEGG gene overlap upset plot"), "upsetplot_enrichKEGG", 12, 8)
  }
}

reactome_result <- tryCatch({
  enrichPathway(
    gene = gene_entrez,
    organism = "human",
    pvalueCutoff = pvalue_cutoff,
    pAdjustMethod = "BH",
    qvalueCutoff = qvalue_cutoff,
    readable = TRUE
  )
}, error = function(e) {
  writeLines(conditionMessage(e), file.path(raw_dir, "ReactomePA_error.txt"))
  NULL
})

if (!is.null(reactome_result)) {
  reactome_table <- save_table(reactome_result, "ReactomePA_enrichPathway_results.csv")
  if (nrow(reactome_table) > 0) {
    safe_plot(dotplot(reactome_result, showCategory = show_category, font.size = 10) + ggtitle("Reactome enrichment - turquoise module"), "dotplot_ReactomePA", 11, 8)
    safe_plot(barplot(reactome_result, showCategory = show_category, font.size = 10) + ggtitle("Reactome enrichment - turquoise module"), "barplot_ReactomePA", 11, 8)
    safe_plot(upsetplot(reactome_result, n = min(15, nrow(reactome_table))) + ggtitle("Reactome gene overlap upset plot"), "upsetplot_ReactomePA", 12, 8)
    reactome_sim <- pairwise_termsim(reactome_result)
    safe_plot(emapplot(reactome_sim, showCategory = min(30, nrow(reactome_table))) + ggtitle("Reactome enrichment map"), "emapplot_ReactomePA", 12, 9)
  }
}

compare_bp <- tryCatch({
  compareCluster(
    geneCluster = list(Turquoise = gene_entrez),
    fun = "enrichGO",
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = pvalue_cutoff,
    qvalueCutoff = qvalue_cutoff,
    readable = TRUE
  )
}, error = function(e) {
  writeLines(conditionMessage(e), file.path(raw_dir, "compareCluster_GO_BP_error.txt"))
  NULL
})

if (!is.null(compare_bp)) {
  compare_table <- save_table(compare_bp, "compareCluster_GO_BP_results.csv")
  if (nrow(compare_table) > 0) {
    safe_plot(dotplot(compare_bp, showCategory = show_category, font.size = 10) + ggtitle("compareCluster GO BP - turquoise module"), "dotplot_compareCluster_GO_BP", 11, 8)
  }
}

all_sources <- list(
  GO_BP = if (!is.null(go_results$BP)) as.data.frame(go_results$BP) else data.frame(),
  GO_MF = if (!is.null(go_results$MF)) as.data.frame(go_results$MF) else data.frame(),
  GO_CC = if (!is.null(go_results$CC)) as.data.frame(go_results$CC) else data.frame(),
  KEGG = if (exists("kegg_readable")) as.data.frame(kegg_readable) else data.frame(),
  Reactome = if (!is.null(reactome_result)) as.data.frame(reactome_result) else data.frame()
)

summary <- purrr::imap_dfr(all_sources, function(tbl, source_name) {
  if (nrow(tbl) == 0) {
    tibble::tibble(source = source_name, significant_terms = 0, top_term = NA_character_, top_adjusted_p = NA_real_)
  } else {
  top <- tbl %>% dplyr::arrange(p.adjust) %>% dplyr::slice(1)
    tibble::tibble(
      source = source_name,
      significant_terms = nrow(tbl),
      top_term = top$Description,
      top_adjusted_p = top$p.adjust
    )
  }
})
write_csv(summary, file.path(processed_dir, "clusterProfiler_source_summary.csv"))

top_overall <- purrr::imap_dfr(all_sources, function(tbl, source_name) {
  if (nrow(tbl) == 0) return(tibble::tibble())
  tbl %>%
    mutate(source = source_name) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = 10)
}) %>%
  dplyr::arrange(p.adjust)
write_csv(top_overall, file.path(processed_dir, "top_enriched_terms_all_sources.csv"))

gsea_note <- paste(
  "GSEA was not run because turquoise_module_genes.csv is an unranked gene list.",
  "clusterProfiler::gseGO and gseKEGG require a named numeric ranked vector,",
  "for example all genes ranked by log2 fold-change."
)
writeLines(gsea_note, file.path(processed_dir, "GSEA_not_run_reason.txt"))

summary_md <- summary %>%
  mutate(
    top_adjusted_p = ifelse(is.na(top_adjusted_p), "", format(top_adjusted_p, scientific = TRUE, digits = 3)),
    top_term = ifelse(is.na(top_term), "", top_term)
  ) %>%
  transmute(line = paste0("| ", source, " | ", significant_terms, " | ", top_term, " | ", top_adjusted_p, " |")) %>%
  pull(line)

top_md <- top_overall %>%
  dplyr::slice_head(n = 15) %>%
  transmute(line = paste0("| ", source, " | ", Description, " | ", Count, " | ", format(p.adjust, scientific = TRUE, digits = 3), " |")) %>%
  pull(line)

report <- paste0(
  "# Turquoise module clusterProfiler analysis\n\n",
  "## Input\n\n",
  "- Input file: inputs/turquoise_module_genes.csv\n",
  "- Cleaned SYMBOL genes: ", nrow(clean_genes), "\n",
  "- Mapped ENTREZ IDs: ", length(gene_entrez), "\n",
  "- Unmapped SYMBOL genes: ", nrow(unmapped), "\n\n",
  "## Methods\n\n",
  "The analysis used clusterProfiler for GO and KEGG over-representation analysis, ReactomePA for Reactome pathway enrichment, and enrichplot for publication-ready visualization. GSEA was not run because the input is an unranked module gene list, not a genome-wide ranked vector.\n\n",
  "## Summary by source\n\n",
  "| Source | Significant terms | Top term | Top adjusted p-value |\n",
  "| --- | ---: | --- | ---: |\n",
  paste(summary_md, collapse = "\n"),
  "\n\n## Top enriched terms\n\n",
  "| Source | Term | Gene count | Adjusted p-value |\n",
  "| --- | --- | ---: | ---: |\n",
  paste(top_md, collapse = "\n"),
  "\n\n## Outputs\n\n",
  "- Results tables are in results/processed/.\n",
  "- PNG plots are in figures/png/.\n",
  "- PDF plots are in figures/pdf/.\n",
  "- Scripts are in scripts/.\n"
)
writeLines(report, file.path(report_dir, "clusterProfiler_turquoise_module_report.md"))

message("clusterProfiler analysis complete.")
message("Clean genes: ", nrow(clean_genes))
message("Mapped Entrez IDs: ", length(gene_entrez))
message("Output folder: ", analysis_dir)
