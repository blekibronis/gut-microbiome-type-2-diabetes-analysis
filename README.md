# Gut Microbiome Type 2 Diabetes Analysis

An R pipeline analyzing gut microbiome differences between Type 2 Diabetes (T2D) and healthy individuals using MGnify 16S rRNA data, alpha/beta diversity, PERMANOVA, and differential abundance (DESeq2, ANCOM-BC2).

## 🏗️ Project Architecture

```text
.
├── README.md                          # Project documentation
├── LICENSE                            # Project license
├── microbiome.Rproj                   # RStudio project file
├── setup_environment.R                # One-time package installation script
├── Bioinformatics Project 1.r         # Main analysis pipeline
├── sampleID.csv                       # User-provided phenotype metadata
├── data_clean/                        # Processed TreeSE / phyloseq objects (.rds)
├── metadata/                          # Cleaned and merged metadata CSVs
├── results/                           # Statistical outputs (.csv, .RData)
├── figures/                           # Generated plots
└── images/                            # Supporting images
```

## 🚀 Getting Started

### 1. Installation (One-time Setup)
Run the setup script to install all required CRAN and Bioconductor packages:
```bash
Rscript setup_environment.R
```

### 2. Run the Analysis
Execute the main pipeline script:
```bash
Rscript "Bioinformatics Project 1.r"
```
