#!/usr/bin/env Rscript
# Differential Proteomics Analysis & DEG-Protein Comparison
# Compare mRNA DEGs with protein expression data
# Requires: DESeq2, limma, ggplot2, igraph

library(ggplot2)
library(tidyverse)
library(gridExtra)

cat("\n", strrep("=", 70), "\n")
cat("Differential Proteomics Analysis Pipeline\n")
cat("mRNA DEGs vs Protein Expression Comparison\n")
cat(strrep("=", 70), "\n\n")

# ============================================================================
# 1. LOAD DATA FUNCTIONS
# ============================================================================

load_degs <- function(file_path) {
  """Load DEG results from your mRNA analysis"""
  cat("[1] Loading DEG data...\n")
  degs <- read.csv(file_path, row.names = 1)
  cat(sprintf("  Loaded %d differentially expressed genes\n", nrow(degs)))
  print(head(degs, 3))
  return(degs)
}

load_proteomics <- function(file_path, data_source = "CPTAC") {
  """Load proteomics data from downloaded source"""
  cat(sprintf("\n[2] Loading %s proteomics data...\n", data_source))
  
  proteins <- read.csv(file_path, row.names = 1)
  cat(sprintf("  Loaded %d proteins x %d samples\n", nrow(proteins), ncol(proteins)))
  
  # Check data structure
  cat(sprintf("  Data type: %s\n", class(proteins[1,1])))
  cat(sprintf("  Value range: [%.2f, %.2f]\n", min(proteins, na.rm=T), max(proteins, na.rm=T)))
  
  return(proteins)
}

# ============================================================================
# 2. GENE-PROTEIN MAPPING
# ============================================================================

map_genes_to_proteins <- function(degs, proteins, gene_col = "gene_name", 
                                   protein_gene_col = NULL) {
  """Map DEG genes to protein data"""
  cat("\n[3] Mapping genes to proteins...\n")
  
  # Extract gene names from DEGs
  if (gene_col %in% colnames(degs)) {
    deg_genes <- toupper(trimws(degs[[gene_col]]))
  } else {
    deg_genes <- toupper(trimws(rownames(degs)))
  }
  
  # Remove empty/NA entries
  deg_genes <- deg_genes[!is.na(deg_genes) & deg_genes != "" & deg_genes != "NA"]
  deg_genes <- unique(deg_genes)
  
  # Extract protein gene names
  if (!is.null(protein_gene_col) && protein_gene_col %in% colnames(proteins)) {
    protein_genes <- toupper(trimws(proteins[[protein_gene_col]]))
  } else {
    # Assume rownames contain gene info
    protein_genes <- toupper(trimws(rownames(proteins)))
  }
  
  protein_genes <- protein_genes[!is.na(protein_genes) & protein_genes != "" & protein_genes != "NA"]
  protein_genes <- unique(protein_genes)
  
  # Find common genes
  common <- intersect(deg_genes, protein_genes)
  
  cat(strrep("-", 70), "\n")
  cat(sprintf("DEG genes: %d\n", length(deg_genes)))
  cat(sprintf("Protein genes: %d\n", length(protein_genes)))
  cat(sprintf("Common (overlap): %d (%.1f%% of DEGs)\n", 
              length(common), (length(common)/length(deg_genes)*100)))
  cat(strrep("-", 70), "\n")
  
  return(list(common = common, 
              deg_genes = deg_genes, 
              protein_genes = protein_genes))
}

# ============================================================================
# 3. CREATE COMPARISON MATRIX
# ============================================================================

create_comparison_matrix <- function(degs, proteins, common_genes, 
                                      gene_col = "gene_name") {
  """Create a comparison matrix of DEGs and proteins"""
  cat("\n[4] Creating comparison matrix...\n")
  
  # Match DEGs
  deg_idx <- which(toupper(trimws(rownames(degs))) %in% toupper(common_genes))
  deg_subset <- degs[deg_idx, ]
  
  # Match proteins
  protein_idx <- which(toupper(trimws(rownames(proteins))) %in% toupper(common_genes))
  protein_subset <- proteins[protein_idx, ]
  
  # Prepare comparison dataframe
  comparison <- data.frame(
    gene_id = rownames(deg_subset),
    log2fc_rna = deg_subset$logFC,
    padj_rna = deg_subset$adj.P.Val,
    stringsAsFactors = FALSE
  )
  
  # Add average protein expression
  protein_means <- rowMeans(protein_subset, na.rm = TRUE)
  comparison$protein_mean <- protein_means[comparison$gene_id]
  
  cat(sprintf("  Comparison matrix: %d genes x %d columns\n", 
              nrow(comparison), ncol(comparison)))
  
  print(head(comparison, 5))
  
  return(comparison)
}

# ============================================================================
# 5. CORRELATION ANALYSIS
# ============================================================================

analyze_correlation <- function(comparison, output_dir = "./results") {
  """Analyze correlation between mRNA and protein changes"""
  cat("\n[6] Correlation analysis...\n")
  
  dir.create(output_dir, showWarnings = FALSE)
  
  # Calculate correlation
  valid <- !is.na(comparison$log2fc_rna) & !is.na(comparison$protein_mean)
  
  if (sum(valid) > 2) {
    corr <- cor.test(
      comparison$log2fc_rna[valid],
      comparison$protein_mean[valid],
      method = "pearson"
    )
    
    cat(sprintf("  mRNA vs Protein Correlation:\n"))
    cat(sprintf("    Pearson r = %.3f\n", corr$estimate))
    cat(sprintf("    p-value = %.2e\n", corr$p.value))
    
    # Create scatter plot
    p <- ggplot(comparison[valid, ], aes(x = log2fc_rna, y = protein_mean)) +
      geom_point(alpha = 0.5, color = "steelblue") +
      geom_smooth(method = "lm", color = "red", alpha = 0.2) +
      theme_minimal() +
      labs(
        title = sprintf("mRNA vs Protein Expression\n(r = %.3f, p = %.2e)", 
                       corr$estimate, corr$p.value),
        x = "mRNA log2 Fold Change",
        y = "Protein Mean Expression"
      )
    
    plot_file <- file.path(output_dir, "correlation_mrna_vs_protein.pdf")
    ggsave(plot_file, p, width = 8, height = 6)
    cat(sprintf("  Plot saved: %s\n", plot_file))
  }
}

# ============================================================================
# 6. GENERATE REPORT
# ============================================================================

generate_report <- function(comparison, common_genes, deg_genes, protein_genes, 
                           output_dir = "./results") {
  """Generate comparison report"""
  cat("\n[7] Generating report...\n")
  
  dir.create(output_dir, showWarnings = FALSE)
  
  report_file <- file.path(output_dir, "DEG_Proteomics_Comparison_Report.txt")
  
  sink(report_file)
  
  cat(strrep("=", 80), "\n")
  cat("mRNA DEGs vs Protein Expression Comparison Report\n")
  cat(strrep("=", 80), "\n\n")
  
  cat(sprintf("Analysis Date: %s\n\n", Sys.time()))
  
  cat("SUMMARY STATISTICS\n")
  cat(strrep("-", 80), "\n")
  cat(sprintf("Total unique DEG genes: %d\n", length(deg_genes)))
  cat(sprintf("Total unique protein genes: %d\n", length(protein_genes)))
  cat(sprintf("Common genes (matched): %d\n", length(common_genes)))
  cat(sprintf("Match percentage: %.1f%% of DEGs found in proteomics\n\n", 
              (length(common_genes)/length(deg_genes)*100)))
  
  cat("TOP 50 DEGs WITH PROTEIN DATA\n")
  cat(strrep("-", 80), "\n")
  
  top <- comparison %>%
    arrange(desc(abs(log2fc_rna))) %>%
    slice(1:50)
  
  for (i in 1:nrow(top)) {
    cat(sprintf("%30s  logFC=%7.2f  padj=%.2e  ProteinExpr=%.2f\n",
                top$gene_id[i],
                top$log2fc_rna[i],
                top$padj_rna[i],
                top$protein_mean[i]))
  }
  
  cat("\n", strrep("=", 80), "\n")
  cat("Report generated by: DEG-Proteomics Analysis Pipeline\n")
  cat(strrep("=", 80), "\n\n")
  
  sink()
  
  cat(sprintf("  Report saved: %s\n", report_file))
}

# ============================================================================
# MAIN EXECUTION
# ============================================================================

main <- function() {
  
  cat("\nUsage Instructions:\n")
  cat(strrep("-", 70), "\n")
  cat("1. Update file paths below with your data locations\n")
  cat("2. Ensure DEG file has columns: gene_name, logFC, adj.P.Val\n")
  cat("3. Protein files should be in standard format (genes x samples)\n")
  cat("4. Run: Rscript differential_proteomics_analysis.R\n\n")
  
  # ========== MODIFY THESE PATHS ==========
  
  deg_file <- "../path/to/your/degs.csv"  # Your mRNA DEG results
  
  proteomics_files <- list(
    CPTAC = "../01_CPTAC/protein_quantification.csv",
    TCGA = "../02_TCGA_GDC/COAD_proteomics.csv",
    TCPA = "../03_TCPA/rppa_data.csv"
  )
  
  output_base <- "./results"
  
  # ========================================
  
  # Check if files exist
  if (!file.exists(deg_file)) {
    cat("ERROR: DEG file not found!\n")
    cat(sprintf("Expected at: %s\n", deg_file))
    cat("Please update the 'deg_file' path in this script.\n")
    return(NULL)
  }
  
  # Load DEG data
  degs <- load_degs(deg_file)
  
  # Process each proteomics source
  for (source in names(proteomics_files)) {
    
    protein_file <- proteomics_files[[source]]
    
    if (!file.exists(protein_file)) {
      cat(sprintf("\nSkipping %s - file not found: %s\n", source, protein_file))
      next
    }
    
    cat(sprintf("\n\n%s\nProcessing %s\n%s\n", 
                strrep("=", 70), source, strrep("=", 70)))
    
    # Create output directory
    output_dir <- file.path(output_base, source)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Load proteomics
    proteins <- load_proteomics(protein_file, source)
    
    # Map genes
    mapping <- map_genes_to_proteins(degs, proteins)
    
    if (length(mapping$common) > 0) {
      # Create comparison
      comparison <- create_comparison_matrix(degs, proteins, mapping$common)
      
      # Analyze correlation
      analyze_correlation(comparison, output_dir)
      
      # Generate report
      generate_report(comparison, mapping$common, mapping$deg_genes, 
                      mapping$protein_genes, output_dir)
    }
  }
  
  cat("\n", strrep("=", 70), "\n")
  cat("Analysis Complete!\n")
  cat("Check 'results/' folder for outputs\n")
  cat(strrep("=", 70), "\n\n")
}

# Run
main()
