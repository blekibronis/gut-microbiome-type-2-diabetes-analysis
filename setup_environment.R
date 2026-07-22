###############################################################################
# ENVIRONMENT SETUP SCRIPT
# Run this script once to install all necessary CRAN and Bioconductor packages.
###############################################################################

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cran_packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "readr",
  "vegan",
  "jsonlite"
)

bioc_packages <- c(
  "MGnifyR",
  "phyloseq",
  "mia",
  "TreeSummarizedExperiment",
  "SummarizedExperiment",
  "S4Vectors",
  "DESeq2",
  "ANCOMBC",
  "microbiome"
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

message("All required packages are installed and ready!")
