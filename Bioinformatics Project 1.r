###############################################################################
# MICROBIOME ANALYSIS PIPELINE: T2D vs. HEALTHY (PRJNA422434 / MGYS00005285)
# 
# Purpose:
#   1. Load phenotype metadata from local sampleID.csv and filter target cohort.
#   2. Query MGnify API for sequence data (taxonomic abundance) & metadata.
#   3. Integrate sequence & phenotype metadata into TreeSE and phyloseq.
#   4. Compute Alpha Diversity (Richness, Shannon, Simpson) and Wilcoxon tests.
#   5. Compute Beta Diversity (Bray-Curtis PCoA, PERMANOVA adonis2, Beta Dispersion).
#   6. Perform Differential Abundance Analysis via DESeq2 & ANCOM-BC2.
#   7. Export all cleaned datasets, statistical reports, and figures.
###############################################################################

# ==============================================================================
# 1. CONFIGURATION & ENVIRONMENT SETUP
# ==============================================================================
SEED <- 123
set.seed(SEED)

STUDY_ID       <- "MGYS00005285"
BIOPROJECT_ID  <- "PRJNA422434"
METADATA_FILE  <- "sampleID.csv"
DISEASE_KEEP   <- c("Healthy", "T2D")
MIN_PREVALENCE <- 0.05 # 5% prevalence threshold for differential abundance

# Create output directories
DIRS <- c("data_raw", "data_clean", "metadata", "results", "figures")
invisible(lapply(DIRS, dir.create, showWarnings = FALSE, recursive = TRUE))

# Load required libraries
suppressPackageStartupMessages({
  library(MGnifyR)
  library(TreeSummarizedExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(mia)
  library(phyloseq)
  library(vegan)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(DESeq2)
  library(ANCOMBC)
})

# Helper function to generate clean taxonomy labels
format_taxa_labels <- function(tax_df) {
  tax_df %>%
    mutate(
      best_taxon_label = case_when(
        !is.na(Species) & Species != "" ~ paste0("s__", Species),
        !is.na(Genus)   & Genus   != "" ~ paste0("g__", Genus),
        !is.na(Family)  & Family  != "" ~ paste0("f__", Family),
        !is.na(Order)   & Order   != "" ~ paste0("o__", Order, "_unclassified"),
        !is.na(Class)   & Class   != "" ~ paste0("c__", Class, "_unclassified"),
        !is.na(Phylum)  & Phylum  != "" ~ paste0("p__", Phylum, "_unclassified"),
        TRUE                            ~ paste0("taxon_", taxon)
      )
    )
}

# ==============================================================================
# 2. LOAD & CLEAN PATIENT METADATA
# ==============================================================================
message("\n--- Step 1: Loading patient phenotype metadata ---")
if (!file.exists(METADATA_FILE)) {
  stop("Metadata file '", METADATA_FILE, "' not found in working directory.")
}

patient_metadata <- read_csv(METADATA_FILE, show_col_types = FALSE) %>%
  filter(BioProject == BIOPROJECT_ID, Disease %in% DISEASE_KEEP) %>%
  mutate(
    disease_group = factor(Disease, levels = DISEASE_KEEP),
    sample_accession_csv = `sample.ID`
  )

write_csv(patient_metadata, "metadata/patient_metadata_PRJNA422434_T2D_Healthy_clean.csv")
message(sprintf("Loaded %d matching samples from patient metadata.", nrow(patient_metadata)))

# ==============================================================================
# 3. QUERY MGNIFY API & INTEGRATE METADATA
# ==============================================================================
message("\n--- Step 2: Fetching sequence data & metadata from MGnify API ---")
mg <- MgnifyClient(useCache = TRUE, cacheDir = ".mgnify_cache")

analyses <- searchAnalysis(mg, type = "studies", accession = STUDY_ID)
mgnify_metadata <- getMetadata(mg, analyses)

write_csv(mgnify_metadata, "metadata/mgnify_metadata_raw.csv")

# Combine patient phenotype metadata with MGnify sample metadata
combined_metadata_clean <- mgnify_metadata %>%
  left_join(patient_metadata, by = c("sample_accession" = "sample_accession_csv")) %>%
  filter(BioProject == BIOPROJECT_ID, Disease %in% DISEASE_KEEP, !is.na(analysis_accession)) %>%
  distinct(analysis_accession, .keep_all = TRUE)

if (nrow(combined_metadata_clean) == 0) {
  stop("No overlapping samples found between MGnify study and sampleID.csv.")
}

write_csv(combined_metadata_clean, "metadata/combined_metadata_PRJNA422434_T2D_Healthy.csv")
selected_analyses <- combined_metadata_clean$analysis_accession

# Automatically remove any 0-byte or corrupted cached BIOM files before downloading
corrupted_cache <- list.files(".mgnify_cache", recursive = TRUE, full.names = TRUE)
corrupted_cache <- corrupted_cache[file.size(corrupted_cache) == 0]
if (length(corrupted_cache) > 0) {
  file.remove(corrupted_cache)
  message(sprintf("Removed %d 0-byte corrupted cache file(s).", length(corrupted_cache)))
}

# Download taxonomic abundance as TreeSummarizedExperiment (TreeSE)
tse <- getResult(
  mg,
  accession  = selected_analyses,
  get.taxa   = "SSU",
  get.func   = FALSE,
  output     = "TreeSE",
  bulk.dl    = FALSE,
  get.tree   = FALSE,
  use.cache  = TRUE
)

# Align metadata and colData
meta_df <- as.data.frame(combined_metadata_clean)
rownames(meta_df) <- meta_df$analysis_accession
common_samples <- intersect(colnames(tse), rownames(meta_df))

tse <- tse[, common_samples]
meta_df <- meta_df[colnames(tse), , drop = FALSE]
colData(tse) <- cbind(colData(tse), S4Vectors::DataFrame(meta_df))

saveRDS(tse, "data_clean/tse_PRJNA422434_T2D_Healthy_with_metadata.rds")

# Convert to Phyloseq and filter low-abundance taxa
ps <- mia::convertToPhyloseq(tse, assay.type = assayNames(tse)[1])
saveRDS(ps, "data_clean/phyloseq_PRJNA422434_T2D_Healthy_with_metadata.rds")

ps_filtered <- prune_taxa(taxa_sums(ps) > 10, ps)
saveRDS(ps_filtered, "data_clean/phyloseq_PRJNA422434_T2D_Healthy_filtered.rds")

message(sprintf("Phyloseq object built: %d samples, %d taxa after low-abundance filtering.", 
                nsamples(ps_filtered), ntaxa(ps_filtered)))

# ==============================================================================
# 4. ALPHA DIVERSITY ANALYSIS
# ==============================================================================
message("\n--- Step 3: Running Alpha Diversity Analysis ---")
alpha_df <- estimate_richness(ps_filtered, measures = c("Observed", "Shannon", "Simpson")) %>%
  mutate(
    sample_key    = rownames(.),
    disease_group = factor(sample_data(ps_filtered)$disease_group, levels = DISEASE_KEEP)
  )

write_csv(alpha_df, "results/alpha_diversity_T2D_vs_Healthy.csv")

# Statistical testing (Wilcoxon Rank-Sum)
alpha_tests <- data.frame(
  metric  = c("Observed", "Shannon", "Simpson"),
  p_value = c(
    wilcox.test(Observed ~ disease_group, data = alpha_df)$p.value,
    wilcox.test(Shannon  ~ disease_group, data = alpha_df)$p.value,
    wilcox.test(Simpson  ~ disease_group, data = alpha_df)$p.value
  )
)
write_csv(alpha_tests, "results/alpha_diversity_wilcox_tests.csv")

# Visualization
p_alpha <- ggplot(alpha_df, aes(x = disease_group, y = Shannon, fill = disease_group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.8) +
  theme_minimal() +
  labs(title = "Alpha Diversity: Shannon Index", x = "Disease Group", y = "Shannon Diversity")

ggsave("figures/alpha_diversity_shannon_T2D_vs_Healthy.png", plot = p_alpha, width = 7, height = 5, dpi = 300)

# ==============================================================================
# 5. BETA DIVERSITY ANALYSIS & PERMANOVA
# ==============================================================================
message("\n--- Step 4: Running Beta Diversity & PERMANOVA ---")
bray_dist <- phyloseq::distance(ps_filtered, method = "bray")

# PCoA Ordination
ord_bray <- ordinate(ps_filtered, method = "PCoA", distance = bray_dist)
p_beta <- plot_ordination(ps_filtered, ord_bray, color = "disease_group") +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(title = "Beta Diversity: Bray-Curtis PCoA", color = "Disease Group")

ggsave("figures/beta_diversity_bray_pcoa_T2D_vs_Healthy.png", plot = p_beta, width = 7, height = 5, dpi = 300)

# PERMANOVA test
permanova_meta <- data.frame(
  disease_group = sample_data(ps_filtered)$disease_group,
  row.names     = sample_names(ps_filtered)
)
sample_order <- labels(bray_dist)
permanova_meta <- permanova_meta[sample_order, , drop = FALSE]

set.seed(SEED)
permanova_result <- vegan::adonis2(
  bray_dist ~ disease_group,
  data         = permanova_meta,
  permutations = 999
)
write.csv(as.data.frame(permanova_result), "results/permanova_bray_T2D_vs_Healthy.csv")

# Beta Dispersion Test
beta_disp <- vegan::betadisper(bray_dist, permanova_meta$disease_group)
png("figures/beta_dispersion_T2D_vs_Healthy.png", width = 800, height = 600)
plot(beta_disp, main = "Beta Dispersion: Healthy vs T2D")
dev.off()

# ==============================================================================
# 6. DIFFERENTIAL ABUNDANCE: DESeq2
# ==============================================================================
message("\n--- Step 5: Running DESeq2 Differential Abundance ---")

# Filter taxa present in >= 5% of samples
min_samples <- ceiling(MIN_PREVALENCE * nsamples(ps_filtered))
otu_mat <- as(otu_table(ps_filtered), "matrix")
if (!taxa_are_rows(ps_filtered)) otu_mat <- t(otu_mat)
taxa_keep <- rowSums(otu_mat > 0) >= min_samples

ps_deseq <- prune_taxa(taxa_keep, ps_filtered)

dds <- phyloseq_to_deseq2(ps_deseq, ~ disease_group)
dds <- DESeq(dds, sfType = "poscounts")
res_deseq <- results(dds, contrast = c("disease_group", "T2D", "Healthy"))

da_deseq <- as.data.frame(res_deseq) %>%
  mutate(taxon_id = rownames(.)) %>%
  left_join(as.data.frame(tax_table(ps_deseq)) %>% mutate(taxon_id = rownames(.)), by = "taxon_id") %>%
  arrange(padj)

write_csv(da_deseq, "results/differential_abundance_DESeq2_T2D_vs_Healthy_poscounts.csv")

top10_deseq <- da_deseq %>% filter(!is.na(padj)) %>% slice_head(n = 10)
write_csv(top10_deseq, "results/top10_biomarker_taxa_DESeq2_poscounts.csv")

# ==============================================================================
# 7. DIFFERENTIAL ABUNDANCE: ANCOM-BC2
# ==============================================================================
message("\n--- Step 6: Running ANCOM-BC2 Differential Abundance ---")

set.seed(SEED)
ancombc2_out <- ancombc2(
  data          = ps_deseq,
  tax_level     = NULL,
  fix_formula   = "disease_group",
  rand_formula  = NULL,
  p_adj_method  = "holm",
  pseudo_sens   = TRUE,
  prv_cut       = MIN_PREVALENCE,
  lib_cut       = 0,
  s0_perc       = 0.05,
  group         = "disease_group",
  struc_zero    = TRUE,
  neg_lb        = TRUE,
  alpha         = 0.05,
  n_cl          = 1,
  verbose       = FALSE,
  global        = FALSE,
  pairwise      = FALSE,
  dunnet        = FALSE,
  trend         = FALSE,
  iter_control  = list(tol = 1e-2, max_iter = 20, verbose = FALSE),
  em_control    = list(tol = 1e-5, max_iter = 100),
  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100)
)

ancom_res <- ancombc2_out$res
write_csv(ancom_res, "results/differential_abundance_ANCOMBC2_raw_T2D_vs_Healthy.csv")

# Extract primary contrast results (T2D vs Healthy reference)
lfc_col  <- grep("^lfc_.*T2D|^lfc_.*disease_group", colnames(ancom_res), value = TRUE)[1]
se_col   <- grep("^se_.*T2D|^se_.*disease_group", colnames(ancom_res), value = TRUE)[1]
w_col    <- grep("^W_.*T2D|^W_.*disease_group", colnames(ancom_res), value = TRUE)[1]
p_col    <- grep("^p_.*T2D|^p_.*disease_group", colnames(ancom_res), value = TRUE)[1]
q_col    <- grep("^q_.*T2D|^q_.*disease_group", colnames(ancom_res), value = TRUE)[1]
diff_col <- grep("^diff_.*T2D|^diff_.*disease_group", colnames(ancom_res), value = TRUE)[1]
rob_col  <- grep("^diff_robust_.*T2D|^diff_robust_.*disease_group", colnames(ancom_res), value = TRUE)[1]

ancom_clean <- ancom_res %>%
  transmute(
    taxon         = taxon,
    lfc           = .data[[lfc_col]],
    se            = if (!is.na(se_col)) .data[[se_col]] else NA_real_,
    W             = if (!is.na(w_col)) .data[[w_col]] else NA_real_,
    p_value       = .data[[p_col]],
    q_value       = .data[[q_col]],
    diff_abundant = .data[[diff_col]],
    diff_robust   = if (!is.na(rob_col)) .data[[rob_col]] else NA
  ) %>%
  mutate(
    direction = case_when(
      lfc > 0 ~ "Higher in T2D",
      lfc < 0 ~ "Higher in Healthy",
      TRUE    ~ "No direction"
    )
  ) %>%
  arrange(q_value)

write_csv(ancom_clean, "results/differential_abundance_ANCOMBC2_clean_T2D_vs_Healthy.csv")

# Attach Taxonomy Annotations
tax_df <- as.data.frame(tax_table(ps_deseq)) %>% mutate(taxon = rownames(.))
ancom_taxonomy <- ancom_clean %>%
  left_join(tax_df, by = "taxon") %>%
  format_taxa_labels()

write_csv(ancom_taxonomy, "results/differential_abundance_ANCOMBC2_with_taxonomy_T2D_vs_Healthy.csv")

# Filter Significant Taxa & Top 10 Biomarkers
ancom_sig <- ancom_taxonomy %>%
  filter(diff_abundant == TRUE | diff_abundant == 1) %>%
  arrange(q_value)

write_csv(ancom_sig, "results/significant_taxa_ANCOMBC2_T2D_vs_Healthy.csv")
write_csv(slice_head(ancom_sig, n = 10), "results/top10_biomarker_taxa_ANCOMBC2.csv")

# Visualization: Top Differentially Abundant Taxa
top_plot_df <- ancom_sig %>%
  slice_head(n = 15) %>%
  mutate(plot_label = make.unique(best_taxon_label))

p_ancom <- ggplot(top_plot_df, aes(x = reorder(plot_label, lfc), y = lfc, fill = direction)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  labs(
    title    = "Top Differentially Abundant Taxa: ANCOM-BC2",
    subtitle = "T2D vs Healthy",
    x        = "Taxon",
    y        = "Log Fold Change (LFC)",
    fill     = "Group Direction"
  )

ggsave("figures/top_differential_abundance_ANCOMBC2_T2D_vs_Healthy.png", plot = p_ancom, width = 8, height = 6, dpi = 300)

# ==============================================================================
# 8. SAVE WORKSPACE & PIPELINE COMPLETE
# ==============================================================================
save.image("results/microbiome_project_workspace.RData")
message("\nPipeline completed successfully! All datasets and figures exported.")
