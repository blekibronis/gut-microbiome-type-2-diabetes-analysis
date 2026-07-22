# Gut Microbiome Type 2 Diabetes Analysis

An R pipeline analyzing gut microbiome differences between Type 2 Diabetes (T2D) and healthy individuals using MGnify 16S rRNA data, alpha/beta diversity, PERMANOVA, and differential abundance (DESeq2, ANCOM-BC2).

## Why this matters?
The gut microbiome has been repeatedly implicated in the pathophysiology of Type 2 Diabetes, particularly through depletion of butyrate-producing bacteria and shifts in a limited set of taxa rather than a wholesale restructuring of the community. This pipeline automates a standard, reproducible microbiome comparison workflow — from MGnify sequence data to statistically tested, taxonomy-annotated biomarker tables.

## 📊 Example Output
 
**Differentially abundant taxa** — log fold change of taxa significantly associated with disease status (ANCOM-BC2, q < 0.05):
 
![Top differentially abundant taxa](/Users/rafaelhutajulu/development/microbiome/images/top_differential_abundance_ANCOMBC2_T2D_vs_Healthy.png)
 
*Butyricicoccus*, a known butyrate producer, is depleted in T2D — consistent with prior literature linking reduced butyrate production to T2D-associated dysbiosis. Enrichment in T2D is concentrated in a small set of taxa (e.g., Christensenellaceae, *Fusobacterium*, *Bilophila*), rather than reflecting a global shift in community composition.
 
**Diversity comparisons**:
 
<table>
<tr>
<td><img src="/Users/rafaelhutajulu/development/microbiome/images/alpha_diversity_shannon_T2D_vs_Healthy.png" alt="Shannon alpha diversity" width="400"/></td>
<td><img src="/Users/rafaelhutajulu/development/microbiome/images/beta_diversity_bray_pcoa_T2D_vs_Healthy.png" alt="Bray-Curtis PCoA" width="400"/></td>
</tr>
</table>
Alpha diversity (Shannon) differs modestly but significantly between groups (Wilcoxon p = 0.036); beta diversity ordination (Bray-Curtis PCoA) shows substantial overlap between T2D and Healthy samples, with PERMANOVA confirming a statistically significant but very small effect (R² = 0.006, p = 0.022).
 
---

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

## 🛠️ Requirements & Installation

### 1. System Requirements
- R version 4.2 or higher
- Internet access for MGnify API queries

### 2. R Package Setup
All required CRAN and Bioconductor packages are installed automatically via the setup script:

```r
source("setup_environment.R")
```

This installs the required packages for data loading, visualization, diversity analysis, and differential abundance testing.

## 🚀 Running the Pipeline

1. Place a sampleID.csv file in the working directory containing at least the columns BioProject, Disease, and sample.ID, matching samples to BioProject PRJNA422434 with disease status Healthy or T2D.
2. Run the pipeline:

```r
source("Bioinformatics Project 1.r")
```

### Configuration
Key parameters are set at the top of the script and can be edited directly:

- STUDY_ID: MGnify study accession to query (default: MGYS00005285)
- BIOPROJECT_ID: BioProject accession used to filter metadata (default: PRJNA422434)
- DISEASE_KEEP: Disease groups to retain (default: c("Healthy", "T2D"))
- MIN_PREVALENCE: Minimum sample prevalence for a taxon to be included in differential abundance testing (default: 0.05)
- SEED: Random seed for reproducibility (default: 123)

The script queries MGnify and caches downloaded results locally in .mgnify_cache/ to avoid re-downloading on subsequent runs.
