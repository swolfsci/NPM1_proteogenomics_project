# =============================================================================
# 00_preprocessing.R
#
# Protein data QC, normalization, and clinical data integration.
#
# Input:  data/protein_intensities_raw.csv
#         data/sample_metadata.csv
#         data/clinical_data.csv
# Output: data/protein_matrix_vsn_imputed.csv  (provided; DreamAI runs on HPC)
#         output/proteomics_qc.RDS
#
# Note: The final analysis matrix (protein_matrix_vsn_imputed.csv) is provided
# in the figshare deposit. This script documents the preprocessing pipeline
# for transparency. DreamAI imputation was performed on an HPC cluster and
# is not re-run here.
# =============================================================================

library(tidyverse)
library(matrixStats)
library(vsn)
library(pcaMethods)
library(patchwork)

source("R/utils.R")

# =============================================================================
# 1. Load raw protein data
# =============================================================================

cat("Loading raw protein data...\n")
prot_data <- read.csv(file.path(data_dir, "protein_intensities_raw.csv"),
                       row.names = 1, check.names = FALSE)

cat(sprintf("Raw data: %d proteins x %d samples\n", nrow(prot_data), ncol(prot_data)))

# =============================================================================
# 2. Quality control: filter samples and proteins with >50% missing values
# =============================================================================

# Remove samples with >50% NA
sample_na_frac <- colSums2(is.na(as.matrix(prot_data))) / nrow(prot_data)
high_na_samples <- names(which(sample_na_frac > 0.5))
cat(sprintf("Removing %d samples with >50%% missing values\n", length(high_na_samples)))
prot_data <- prot_data[, !colnames(prot_data) %in% high_na_samples]

# Remove proteins with >50% NA
protein_na_frac <- rowSums2(is.na(as.matrix(prot_data))) / ncol(prot_data)
high_na_proteins <- names(which(protein_na_frac > 0.5))
cat(sprintf("Removing %d proteins with >50%% missing values\n", length(high_na_proteins)))
prot_data <- prot_data[!rownames(prot_data) %in% high_na_proteins, ]

cat(sprintf("After QC: %d proteins x %d samples\n", nrow(prot_data), ncol(prot_data)))

# =============================================================================
# 3. VSN normalization + median centering
# =============================================================================

cat("Performing VSN normalization...\n")
prot_data_vsn <- justvsn(as.matrix(prot_data))

# Median centering to correct for residual loading differences
vsn_col_medians <- colMedians(prot_data_vsn, na.rm = TRUE)
prot_data_vsn_centered <- sweep(prot_data_vsn, 2, vsn_col_medians, "-")

cat("VSN + median centering complete.\n")

# =============================================================================
# 4. DreamAI imputation (performed on HPC cluster)
# =============================================================================

# The VSN-normalized, median-centered matrix was exported and imputed using
# DreamAI (Ma et al., bioRxiv 2021) with default settings on an HPC cluster.
# The imputed matrix is provided as: data/protein_matrix_vsn_imputed.csv
#
# To reproduce the imputation:
#   library(DreamAI)
#   result <- DreamAI(prot_data_vsn_centered, k = 10, maxiter_MF = 10,
#                     ntree = 100, maxnodes = NULL, maxiter_ADMIN = 30,
#                     tol = 10^(-2), gamma_ADMIN = 0, gamma = 50,
#                     CV = FALSE, fillmethod = "row_mean", maxiter_RegImpute = 10,
#                     conv_nrmse = 1e-06, iter_SpectroFM = 40, method = c(
#                       "KNN", "MissForest", "ADMIN", "Brinn", "SpectroFM",
#                       "RegImpute"), out = "Ensemble")
#   prot_data_imputed <- result$Ensemble

cat("Loading pre-computed imputed matrix...\n")
prot_data_imputed <- as.matrix(
  read.csv(file.path(data_dir, "protein_matrix_vsn_imputed.csv"),
           row.names = 1, check.names = FALSE)
)

# Remove internal standard (IS) samples, keep only patient samples
is_patient <- grepl("^bioid_", colnames(prot_data_imputed))
prot_data_aml <- prot_data_imputed[, is_patient]

cat(sprintf("Final analysis matrix: %d proteins x %d patients\n",
            nrow(prot_data_aml), ncol(prot_data_aml)))

# =============================================================================
# 5. Clinical data integration
# =============================================================================

cat("Loading clinical data...\n")
clinical_data <- read.csv(file.path(data_dir, "clinical_data.csv"))

clinical_data <- clinical_data %>%
  mutate(
    bio_id_merge = paste0("bioid_", bio_id),
    efs_months = efs_days %/% 30.4375,
    os_months = os_days %/% 30.4375,
    LDH_diag_log10 = log10(LDH_diag),
    wbc_log10 = log10(wbc),
    cohesin = as.numeric(RAD21 == 1 | SMC1A == 1 | SMC3 == 1 | STAG2 == 1),
    ras = as.numeric(NRAS == 1 | KRAS == 1),
    DNMT3A_R882 = as.numeric(DNMT3A_Hotspot == "R882" & !is.na(DNMT3A_Hotspot)),
    age_cat = case_when(
      age < 50 ~ "<50",
      age >= 50 & age < 70 ~ ">=50 & <70",
      age >= 70 ~ ">=70"
    )
  ) %>%
  filter(bio_id_merge %in% colnames(prot_data_aml))

cat(sprintf("Clinical data for %d patients\n", nrow(clinical_data)))

# =============================================================================
# 6. QC visualization
# =============================================================================

# PCA before and after normalization for QC
pca_vsn <- prcomp(t(prot_data_aml), center = TRUE, scale. = FALSE)
var_explained <- summary(pca_vsn)$importance[2, 1:5] * 100

p_qc <- tibble(
  bio_id_merge = rownames(pca_vsn$x),
  PC1 = pca_vsn$x[, 1],
  PC2 = pca_vsn$x[, 2]
) %>%
  ggplot(aes(PC1, PC2)) +
  geom_point(size = 1.5, alpha = 0.7) +
  labs(
    x = sprintf("PC1 (%.1f%%)", var_explained[1]),
    y = sprintf("PC2 (%.1f%%)", var_explained[2]),
    title = "PCA of normalized, imputed proteome"
  ) +
  cowplot::theme_cowplot()

ggsave(file.path(fig_dir, "qc_pca.pdf"), p_qc, width = 6, height = 5)

# =============================================================================
# 7. Save intermediate objects for downstream scripts
# =============================================================================

saveRDS(list(
  vsn_matrix = prot_data_aml,
  clinical = clinical_data
), file.path(output_dir, "preprocessed_data.RDS"))

cat("Preprocessing complete. Saved: output/preprocessed_data.RDS\n")
