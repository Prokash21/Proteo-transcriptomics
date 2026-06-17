#!/usr/bin/env python3
"""
Multi-source Proteomics Data Downloader
Downloads proteomics data for colon and rectum cancer from multiple sources
Organizes data in separate folders for each source

Author: Analysis Pipeline
Date: 2026
"""

import os
import sys
import time
import logging
from pathlib import Path
from datetime import datetime

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('proteomics_download.log', encoding='utf-8'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

# Define base directories
BASE_DIR = Path(__file__).parent.parent
SOURCES = {
    'CPTAC': BASE_DIR / '01_CPTAC',
    'TCGA_GDC': BASE_DIR / '02_TCGA_GDC',
    'TCPA': BASE_DIR / '03_TCPA',
    'ProteomicsDB': BASE_DIR / '04_ProteomicsDB'
}

CANCER_TYPES = {
    'COAD': 'Colon Adenocarcinoma',
    'READ': 'Rectal Adenocarcinoma'
}


class CPTAC_Downloader:
    """Download data from CPTAC (Clinical Proteomic Tumor Analysis Consortium)"""
    
    def __init__(self, output_dir):
        self.output_dir = output_dir
        self.base_url = "https://proteomics.cancer.gov"
        
    def download(self):
        logger.info("=" * 60)
        logger.info("CPTAC Data Download")
        logger.info("=" * 60)
        logger.info("CPTAC offers the most comprehensive proteomics data.")
        logger.info("Manual download recommended due to access restrictions.")
        logger.info("")
        logger.info("Steps to download CPTAC data:")
        logger.info("1. Visit: https://proteomics.cancer.gov/")
        logger.info("2. Go to 'Browse Data' → Select colorectal cancer datasets")
        logger.info("3. Look for datasets like:")
        logger.info("   - CPTAC Colorectal Cancer (COAD)")
        logger.info("   - Download 'Protein quantification' matrices")
        logger.info("4. Save files to: " + str(self.output_dir))
        logger.info("")
        logger.info("Expected file formats: .txt, .csv (log2-normalized or raw)")
        logger.info("")


class TCGA_GDC_Downloader:
    """Download data from GDC (Genomic Data Commons) TCGA"""
    
    def __init__(self, output_dir):
        self.output_dir = output_dir
        self.api_url = "https://api.gdc.cancer.gov"
        
    def create_query_script(self):
        """Create an R script for automated GDC download using TCGAbiolinks"""
        script_content = '''#!/usr/bin/env Rscript
# Download TCGA Proteomics Data using TCGAbiolinks
# For Colon (COAD) and Rectum (READ) cancer

# Install packages if needed
# install.packages("TCGAbiolinks")
# install.packages("SummarizedExperiment")

library(TCGAbiolinks)
library(SummarizedExperiment)

# Download Protein Expression for COAD (Colon)
cat("\\nDownloading TCGA-COAD Protein Expression...\\n")

query_COAD <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Protein Expression",
  data.type = "Protein Expression Quantification",
  workflow.type = "CDAP"
)

# GDCdownload(query_COAD, directory = "./02_TCGA_GDC/COAD")
# data_COAD <- GDCprepare(query_COAD, directory = "./02_TCGA_GDC/COAD")
# save(data_COAD, file = "./02_TCGA_GDC/COAD/TCGA_COAD_Proteomics.RData")

# Download Protein Expression for READ (Rectum)
cat("\\nDownloading TCGA-READ Protein Expression...\\n")

query_READ <- GDCquery(
  project = "TCGA-READ",
  data.category = "Protein Expression",
  data.type = "Protein Expression Quantification",
  workflow.type = "CDAP"
)

# GDCdownload(query_READ, directory = "./02_TCGA_GDC/READ")
# data_READ <- GDCprepare(query_READ, directory = "./02_TCGA_GDC/READ")
# save(data_READ, file = "./02_TCGA_GDC/READ/TCGA_READ_Proteomics.RData")

cat("\\nDownload script created. Uncomment GDCdownload() lines to execute.\\n")
cat("Note: First downloads require valid GDC authentication credentials.\\n")
'''
        
        script_path = self.output_dir / "download_tcga_proteomics.R"
        with open(script_path, 'w') as f:
            f.write(script_content)
        
        return script_path
    
    def download(self):
        logger.info("=" * 60)
        logger.info("TCGA/GDC Data Download")
        logger.info("=" * 60)
        logger.info("Creating R script for automated download via TCGAbiolinks...")
        
        script_path = self.create_query_script()
        logger.info(f"R script created: {script_path}")
        logger.info("")
        logger.info("To download TCGA proteomics data:")
        logger.info("1. Install R packages: install.packages('TCGAbiolinks')")
        logger.info("2. Run: Rscript 02_TCGA_GDC/download_tcga_proteomics.R")
        logger.info("3. Or manually via: https://portal.gdc.cancer.gov/")
        logger.info("   - Project: TCGA-COAD or TCGA-READ")
        logger.info("   - Data Category: Protein Expression")
        logger.info("")


class TCPA_Downloader:
    """Download data from TCPA (The Cancer Protein Atlas)"""
    
    def __init__(self, output_dir):
        self.output_dir = output_dir
        self.base_url = "https://tcpaportal.org"
        
    def create_query_script(self):
        """Create Python script to query TCPA API"""
        script_content = '''#!/usr/bin/env python3
"""
Download TCPA (The Cancer Protein Atlas) data
Requires: pandas, requests
"""

import requests
import pandas as pd
from pathlib import Path

BASE_DIR = Path(__file__).parent.parent / "03_TCPA"

def download_tcpa_data():
    """Query TCPA API for colorectal cancer protein expression"""
    
    # TCPA provides pre-compiled data for different cancer types
    # Manual download: https://tcpaportal.org/
    
    print("TCPA Data Download Instructions")
    print("=" * 50)
    print("TCPA has curated protein expression data from TCGA")
    print()
    print("Steps:")
    print("1. Visit: https://tcpaportal.org/")
    print("2. Select Cancer Type: Colorectal Cancer (COAD/READ)")
    print("3. Download protein expression matrix")
    print("4. Save to: " + str(BASE_DIR))
    print()
    print("Note: TCPA data is already normalized and processed")
    print("Includes RPPA (Reverse Phase Protein Array) data")
    print()

if __name__ == "__main__":
    download_tcpa_data()
'''
        
        script_path = self.output_dir / "download_tcpa_data.py"
        with open(script_path, 'w') as f:
            f.write(script_content)
        
        return script_path
    
    def download(self):
        logger.info("=" * 60)
        logger.info("TCPA Data Download")
        logger.info("=" * 60)
        logger.info("Creating query script for TCPA...")
        
        script_path = self.create_query_script()
        logger.info(f"Python script created: {script_path}")
        logger.info("")
        logger.info("TCPA Manual Download Steps:")
        logger.info("1. Visit: https://tcpaportal.org/")
        logger.info("2. Select cancer type: Colorectal Cancer")
        logger.info("3. Download protein expression data")
        logger.info("4. Save to: 03_TCPA folder")
        logger.info("")
        logger.info("TCPA Strengths:")
        logger.info("  - Pre-normalized RPPA data")
        logger.info("  - Ready for comparative analysis")
        logger.info("  - Includes antibody information")
        logger.info("")


class ProteomicsDB_Downloader:
    """Download data from ProteomicsDB"""
    
    def __init__(self, output_dir):
        self.output_dir = output_dir
        self.base_url = "https://www.proteomicsdb.org"
        
    def create_query_script(self):
        """Create Python script to query ProteomicsDB API"""
        script_content = '''#!/usr/bin/env python3
"""
Download ProteomicsDB data
Requires: requests, pandas
"""

import requests
import pandas as pd
from pathlib import Path
import json

BASE_DIR = Path(__file__).parent.parent / "04_ProteomicsDB"

def query_proteomicsdb():
    """Query ProteomicsDB API for protein expression"""
    
    # ProteomicsDB API endpoints
    api_base = "https://www.proteomicsdb.org/api/v1"
    
    print("ProteomicsDB Data Access Instructions")
    print("=" * 50)
    print()
    print("Option 1: Manual Web Download")
    print("-" * 50)
    print("1. Visit: https://www.proteomicsdb.org/")
    print("2. Search for 'colorectal cancer'")
    print("3. Download tissue-specific protein expression")
    print()
    
    print("Option 2: API Access")
    print("-" * 50)
    print("Example to get proteins:")
    print()
    
    # Example: Query for specific proteins
    example_url = api_base + "/protein/search/name/TP53"
    print(f"  URL: {example_url}")
    print()
    
    print("Save results to: " + str(BASE_DIR))

if __name__ == "__main__":
    query_proteomicsdb()
'''
        
        script_path = self.output_dir / "download_proteomicsdb_data.py"
        with open(script_path, 'w') as f:
            f.write(script_content)
        
        return script_path
    
    def download(self):
        logger.info("=" * 60)
        logger.info("ProteomicsDB Data Download")
        logger.info("=" * 60)
        logger.info("Creating query script for ProteomicsDB...")
        
        script_path = self.create_query_script()
        logger.info(f"Python script created: {script_path}")
        logger.info("")
        logger.info("ProteomicsDB Access Options:")
        logger.info("1. Web Portal: https://www.proteomicsdb.org/")
        logger.info("2. REST API: Programmatic data access")
        logger.info("3. Download tissue-specific proteome data")
        logger.info("4. Save to: 04_ProteomicsDB folder")
        logger.info("")


def create_summary_document():
    """Create a comprehensive guide document"""
    
    guide_content = '''# Proteomics Data Download Guide
## For Differential Proteomics Analysis: Colon & Rectum Cancer

Generated: {timestamp}

---

## Data Sources Overview

### 1. CPTAC (01_CPTAC/)
**Best for:** Comprehensive mass spectrometry proteomics

**URL:** https://proteomics.cancer.gov/

**What to download:**
- Protein quantification matrices
- Label-free MS/MS data
- Log2-normalized or raw counts

**Process:**
1. Browse datasets → Colorectal Cancer
2. Download "Protein quantification" files
3. Save to 01_CPTAC/

**File format:** .txt, .csv (tab/comma-delimited)

---

### 2. TCGA/GDC (02_TCGA_GDC/)
**Best for:** Well-annotated clinical data with proteomics

**URL:** https://portal.gdc.cancer.gov/

**Projects to search:**
- TCGA-COAD (Colon Adenocarcinoma)
- TCGA-READ (Rectal Adenocarcinoma)

**Filter by:**
- Data Category: Protein Expression
- Data Type: Protein Expression Quantification

**Access methods:**
- Web portal (manual)
- `gdc-client` (command-line)
- R: `TCGAbiolinks` package
- Python: `GDC API`

**Quick R script:** See `download_tcga_proteomics.R`

---

### 3. TCPA (03_TCPA/)
**Best for:** RPPA data (curated from TCGA)

**URL:** https://tcpaportal.org/

**Advantages:**
- Pre-normalized RPPA (Reverse Phase Protein Array)
- Quality controlled
- Ready for analysis
- Includes antibody information

**Download:**
1. Select Cancer Type: Colorectal
2. Download protein expression matrix
3. Save to 03_TCPA/

---

### 4. ProteomicsDB (04_ProteomicsDB/)
**Best for:** Tissue-specific proteome reference

**URL:** https://www.proteomicsdb.org/

**Features:**
- Tissue-specific protein abundance
- MS/MS spectral data
- Gene/Protein mapping

**Search:** 
- Filter by tissue/cancer type
- Download expression profiles

---

## Recommended Download Strategy

### Priority Order:
1. **CPTAC** - Most comprehensive (mass spec)
2. **TCGA/GDC** - Most samples, clinical data
3. **TCPA** - Pre-processed, normalized
4. **ProteomicsDB** - Reference/validation data

---

## Next Steps After Download

### 1. Data Exploration (Python/R)
```python
import pandas as pd

# Load proteomics data
proteins = pd.read_csv('protein_expression.csv', index_col=0)
print(proteins.shape)  # (proteins, samples)
print(proteins.head())
```

### 2. Normalize & Process
- Check for missing values
- Normalize across samples
- Log2 transformation if needed
- Quality control metrics

### 3. Match with DEGs
- Map protein IDs to gene names
- Compare with your mRNA DEG list
- Create overlap matrix
- Prepare for correlation analysis

### 4. Differential Proteomics Analysis
- DEP or limma package (R)
- scipy.stats or statsmodels (Python)
- Design matrix: case vs control
- Statistical testing: t-test or linear model
- Multiple testing correction (BH, Benjamini-Hochberg)

---

## File Organization

```
Proteomics/
├── 01_CPTAC/
│   ├── protein_quantification.csv
│   ├── sample_metadata.csv
│   └── README.txt
│
├── 02_TCGA_GDC/
│   ├── COAD/
│   │   └── TCGA_COAD_Proteomics.RData
│   ├── READ/
│   │   └── TCGA_READ_Proteomics.RData
│   └── download_tcga_proteomics.R
│
├── 03_TCPA/
│   ├── COAD_RPPA.csv
│   ├── READ_RPPA.csv
│   └── antibody_info.txt
│
├── 04_ProteomicsDB/
│   ├── colorectal_proteome.csv
│   └── protein_metadata.txt
│
└── scripts/
    ├── download_proteomics_data.py (this file)
    ├── download_tcga_proteomics.R
    ├── download_tcpa_data.py
    └── download_proteomicsdb_data.py
```

---

## Key Considerations

### Data Formats
- CSV/TSV: Rows = proteins, Columns = samples
- RData (R): R objects, requires R to process
- HDF5: Large datasets, hierarchical

### ID Mapping
- Gene ID vs Protein ID vs UniProt ID
- Use tools: biomaRt (R), mygene (Python)
- Create mapping file for cross-reference

### Normalization
- CPTAC: Often log2-normalized
- TCGA: May need normalization
- TCPA: Pre-normalized
- Always check documentation

### Quality Control
- Sample coverage (% proteins detected)
- Missing values handling
- Outlier detection
- PCA/heatmap visualization

---

## Useful R Packages

```R
# Data retrieval
install.packages("TCGAbiolinks")
install.packages("SummarizedExperiment")

# Processing
install.packages("limma")
install.packages("DEP")
install.packages("tidyverse")

# Visualization
install.packages("ggplot2")
install.packages("pheatmap")
install.packages("ComplexHeatmap")

# Analysis
install.packages("igraph")
install.packages("STRINGdb")
```

---

## Useful Python Packages

```bash
pip install pandas numpy scipy scikit-learn
pip install matplotlib seaborn plotly
pip install biomart mygene
pip install statsmodels scikit-postnos
pip install requests
```

---

## References

1. CPTAC: https://proteomics.cancer.gov/
2. TCGA/GDC: https://portal.gdc.cancer.gov/
3. TCPA: https://tcpaportal.org/
4. ProteomicsDB: https://www.proteomicsdb.org/
5. TCGAbiolinks: https://bioconductor.org/packages/TCGAbiolinks/
6. DEP paper: Proteomics. 2017;17(3-4):1600025

---

**Last Updated:** {timestamp}
**Status:** Ready for data download
'''
    
    guide_path = BASE_DIR / '00_DOWNLOAD_GUIDE.md'
    with open(guide_path, 'w', encoding='utf-8') as f:
        f.write(guide_content.format(timestamp=datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    
    return guide_path


def main():
    """Main execution function"""
    
    logger.info("\n" + "=" * 70)
    logger.info("PROTEOMICS DATA DOWNLOAD ORCHESTRATOR")
    logger.info("Multi-source Colorectal Cancer Proteomics Data Retrieval")
    logger.info("=" * 70 + "\n")
    
    try:
        # Verify output directories
        for source, path in SOURCES.items():
            if not path.exists():
                path.mkdir(parents=True, exist_ok=True)
                logger.info(f"Created directory: {path}")
        
        # Initialize downloaders
        logger.info("\nInitializing data downloaders...\n")
        
        downloaders = {
            'CPTAC': CPTAC_Downloader(SOURCES['CPTAC']),
            'TCGA_GDC': TCGA_GDC_Downloader(SOURCES['TCGA_GDC']),
            'TCPA': TCPA_Downloader(SOURCES['TCPA']),
            'ProteomicsDB': ProteomicsDB_Downloader(SOURCES['ProteomicsDB']),
        }
        
        # Execute downloads
        for name, downloader in downloaders.items():
            downloader.download()
            time.sleep(1)
        
        # Create comprehensive guide
        logger.info("=" * 60)
        logger.info("Creating comprehensive download guide...")
        guide_path = create_summary_document()
        logger.info(f"Guide created: {guide_path}")
        logger.info("=" * 60)
        
        # Final summary
        logger.info("\n" + "=" * 70)
        logger.info("DOWNLOAD PREPARATION COMPLETE!")
        logger.info("=" * 70)
        logger.info("\nNext Steps:")
        logger.info("1. Read: 00_DOWNLOAD_GUIDE.md (comprehensive instructions)")
        logger.info("2. Follow source-specific instructions above")
        logger.info("3. Download data to respective folders")
        logger.info("4. Run preprocessing scripts once data is downloaded")
        logger.info("\nData will be organized as:")
        for source, path in SOURCES.items():
            logger.info(f"  {source:20} → {path}")
        logger.info("\n" + "=" * 70 + "\n")
        
    except Exception as e:
        logger.error(f"Error during execution: {str(e)}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
