#!/usr/bin/env python3
"""
Comparison Script: mRNA DEGs vs Protein Expression
Matches DEGs with proteomics data and identifies common proteins
"""

import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict

class DEG_Proteomics_Comparator:
    """Compare mRNA DEGs with protein expression data"""
    
    def __init__(self, deg_file, proteomics_dir):
        """
        Args:
            deg_file: Path to your DEG results (CSV with gene names/IDs)
            proteomics_dir: Path to proteomics data directory
        """
        self.deg_file = Path(deg_file)
        self.proteomics_dir = Path(proteomics_dir)
        self.results = {}
        
    def load_deg_data(self):
        """Load your DEG data"""
        try:
            self.degs = pd.read_csv(self.deg_file)
            print(f"Loaded DEGs: {self.degs.shape[0]} genes")
            print(f"Columns: {list(self.degs.columns)}")
            return self.degs
        except Exception as e:
            print(f"Error loading DEG file: {e}")
            return None
    
    def load_proteomics_data(self, source='CPTAC'):
        """Load proteomics data from a specific source"""
        
        source_paths = {
            'CPTAC': self.proteomics_dir / '01_CPTAC' / 'protein_quantification.csv',
            'TCGA': self.proteomics_dir / '02_TCGA_GDC' / 'protein_data.csv',
            'TCPA': self.proteomics_dir / '03_TCPA' / 'rppa_data.csv',
            'ProteomicsDB': self.proteomics_dir / '04_ProteomicsDB' / 'proteome_data.csv',
        }
        
        protein_file = source_paths.get(source)
        
        if not protein_file or not protein_file.exists():
            print(f"Protein file not found: {protein_file}")
            return None
        
        try:
            proteins = pd.read_csv(protein_file, index_col=0)
            print(f"Loaded {source} proteins: {proteins.shape[0]} proteins")
            return proteins
        except Exception as e:
            print(f"Error loading {source} proteomics: {e}")
            return None
    
    def map_genes_to_proteins(self, degs, proteins, gene_col='gene_name', protein_col=None):
        """
        Match genes from DEGs to proteins
        
        Args:
            degs: DEG dataframe
            proteins: Protein expression dataframe
            gene_col: Column name with gene identifiers
            protein_col: Column with protein IDs (auto-detect if None)
        """
        
        # Get gene names from DEGs, clean up
        deg_genes = set(degs[gene_col].str.strip().str.upper().unique())
        deg_genes = {g for g in deg_genes if g and g != 'NA' and g != 'NAN'}
        
        # Get protein gene names (index or column)
        if protein_col:
            protein_genes = set(proteins[protein_col].str.strip().str.upper().unique())
        else:
            # Assume index contains gene names
            protein_genes = set(proteins.index.str.strip().str.upper().unique())
        
        protein_genes = {p for p in protein_genes if p and p != 'NA' and p != 'NAN'}
        
        # Find intersection
        common_genes = deg_genes.intersection(protein_genes)
        
        print(f"\n{'='*60}")
        print(f"Gene/Protein Matching Results:")
        print(f"{'='*60}")
        print(f"Total DEG genes: {len(deg_genes)}")
        print(f"Total proteins with gene info: {len(protein_genes)}")
        print(f"Common (DEG & Proteomics): {len(common_genes)}")
        print(f"Match percentage: {(len(common_genes)/len(deg_genes)*100):.1f}%")
        print(f"{'='*60}\n")
        
        return common_genes, deg_genes, protein_genes
    
    def create_comparison_matrix(self, degs, proteins, common_genes, gene_col='gene_name'):
        """Create a comparison matrix of DEGs and protein expression"""
        
        # Filter to common genes
        deg_subset = degs[degs[gene_col].str.upper().isin(common_genes)].copy()
        
        # Create results dataframe
        comparison = pd.DataFrame({
            'gene_id': deg_subset[gene_col],
            'log2fc_rna': deg_subset.get('log2FoldChange', deg_subset.get('logFC', np.nan)),
            'padj_rna': deg_subset.get('padj', deg_subset.get('adj.P.Val', np.nan)),
        })
        
        print(f"Comparison matrix shape: {comparison.shape}")
        print(f"Columns: {list(comparison.columns)}")
        
        return comparison
    
    def generate_report(self, comparison, output_dir='./'):
        """Generate summary report"""
        
        report_path = Path(output_dir) / 'DEG_Proteomics_Comparison_Report.txt'
        
        with open(report_path, 'w') as f:
            f.write("="*70 + "\n")
            f.write("mRNA DEGs vs Protein Expression Comparison Report\n")
            f.write("="*70 + "\n\n")
            
            f.write(f"Total genes in comparison: {len(comparison)}\n")
            f.write(f"Genes with DEG data: {comparison['log2fc_rna'].notna().sum()}\n")
            f.write(f"UniProt ID ratio: {(comparison['protein_id'].notna().sum() / len(comparison) * 100):.1f}%\n\n")
            
            f.write("Top 50 genes with mRNA changes:\n")
            f.write("-"*70 + "\n")
            
            top_genes = comparison.nlargest(50, 'log2fc_rna')
            for idx, row in top_genes.iterrows():
                f.write(f"{row['gene_id']:20} log2FC={row['log2fc_rna']:7.2f} padj={row['padj_rna']:.2e}\n")
            
        print(f"Report saved: {report_path}")
        return report_path


def main():
    """Example usage"""
    
    print("\n" + "="*70)
    print("mRNA DEGs vs Proteomics Comparison Tool")
    print("="*70 + "\n")
    
    # Paths (modify for your data)
    dag_file = Path("../your_degs.csv")  # Your DEG results
    proteomics_dir = Path("../")
    
    comparator = DEG_Proteomics_Comparator(dag_file, proteomics_dir)
    
    # Load data
    degs = comparator.load_deg_data()
    if degs is None:
        print("ERROR: Could not load DEG file. Update path in script.")
        print(f"Expected at: {dag_file}")
        return
    
    # Try loading from each source
    for source in ['CPTAC', 'TCGA', 'TCPA', 'ProteomicsDB']:
        print(f"\nProcessing {source} data...")
        proteins = comparator.load_proteomics_data(source)
        
        if proteins is not None:
            # Perform comparison
            common, deg_set, protein_set = comparator.map_genes_to_proteins(
                degs, proteins,
                gene_col='gene_name'  # Adjust to match your DEG file
            )
            
            if common:
                comparison = comparator.create_comparison_matrix(degs, proteins, common)
                comparator.generate_report(comparison, output_dir=f'./results_{source}')


if __name__ == "__main__":
    print("\n" + "="*70)
    print("Usage Instructions:")
    print("="*70)
    print("""
1. Update the 'dag_file' variable to point to your DEG results CSV
2. Ensure proteomics data is downloaded to appropriate folders
3. Run this script to compare DEGs with protein expression
4. Results will be saved to results_[SOURCE]/ folders

Your DEG file should have columns like:
  - gene_name (or gene_id)
  - log2FoldChange (or logFC)
  - padj (adjusted p-value)

Example:
  python compare_degs_proteomics.py
""")
    print("="*70 + "\n")
    
    main()
