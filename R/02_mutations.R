# =============================================================================
# 02_mutations.R
#
# Mutational landscape and genotype-phenotype associations.
# Generates: Figures 2 and S2 (remaining panels)
#
# Input:  output/preprocessed_data.RDS
#         output/differentiation_data.RDS
# Output: output/figures/Fig2_*.pdf
#         output/figures/FigS2_*.pdf
# =============================================================================

library(tidyverse)
library(ComplexHeatmap)
library(survival)
library(survminer)

source("R/utils.R")

# =============================================================================
# 1. Load data
# =============================================================================

cat("Loading data...\n")
preprocessed <- readRDS(file.path(output_dir, "preprocessed_data.RDS"))
vsn <- preprocessed$vsn_matrix
clinical <- preprocessed$clinical

diff_data <- readRDS(file.path(output_dir, "differentiation_data.RDS"))
cluster_mapping <- diff_data$cluster_mapping
dm_coords <- diff_data$dm_coords
zeng_gsva_long <- diff_data$zeng_gsva_long
diff_scores_wide <- diff_data$diff_scores_wide

# =============================================================================
# 2. OncoPrint (Fig 2A)
# =============================================================================

# Prepare mutation matrix ordered by cluster
ngs_cols <- clinical %>%
  dplyr::select(bio_id_merge, FLT3_ITD_PCR, NRAS, KRAS, DNMT3A_R882,
                PTPN11, RAD21, SMC1A, SMC3, STAG2, IDH1, IDH2, TET2,
                DNMT3A, WT1, RUNX1, ASXL1, ras, cohesin) %>%
  distinct(bio_id_merge, .keep_all = TRUE)

ngs_mat <- cluster_mapping %>%
  dplyr::select(bio_id_merge, cluster) %>%
  left_join(ngs_cols, by = "bio_id_merge") %>%
  arrange(match(cluster, c("Immature_like", "GMP_like", "Commited_like", "Intermediate"))) %>%
  dplyr::select(-cluster) %>%
  column_to_rownames("bio_id_merge") %>%
  as.matrix()

# Filter to mutations with >=3% frequency
mut_freq <- colSums(ngs_mat, na.rm = TRUE) / nrow(ngs_mat) * 100
ngs_mat_filt <- ngs_mat[, mut_freq >= 3]

# Convert to oncoprint format
ngs_mat_onco <- apply(ngs_mat_filt, 2, function(x) ifelse(x == 0 | is.na(x), "", "mutation"))
ngs_mat_onco <- t(ngs_mat_onco)

alter_fun <- list(
  background = function(x, y, w, h) grid::grid.rect(x, y, w, h,
    gp = grid::gpar(fill = "#F0F0F0", col = NA)),
  mutation = function(x, y, w, h) grid::grid.rect(x, y, w * 0.9, h * 0.9,
    gp = grid::gpar(fill = "#e63946", col = NA))
)

ordered_clusters <- cluster_mapping %>%
  arrange(match(cluster, c("Immature_like", "GMP_like", "Commited_like", "Intermediate")))

oncoprint_top_anno <- HeatmapAnnotation(
  df = data.frame(cluster = ordered_clusters$cluster),
  col = list(cluster = cluster_colors)
)

pdf(file.path(fig_dir, "Fig2A_oncoprint.pdf"), width = 12, height = 6)
oncoPrint(
  ngs_mat_onco[, ordered_clusters$bio_id_merge],
  alter_fun = alter_fun,
  col = c("mutation" = "#e63946"),
  show_column_names = FALSE,
  show_row_names = TRUE,
  remove_empty_columns = TRUE,
  remove_empty_rows = TRUE,
  column_title = "Mutation Landscape",
  column_split = ordered_clusters$cluster,
  top_annotation = oncoprint_top_anno
)
dev.off()

# =============================================================================
# 3. Mutation frequency by cluster (Fig 2B)
# =============================================================================

# Fisher's exact test for mutation enrichment per cluster
sig_diff_muts <- c("FLT3_ITD_PCR", "ras", "DNMT3A_R882", "PTPN11", "RAD21")

cluster_sizes <- cluster_mapping %>%
  summarize(n = n(), .by = cluster)

p_mut_freq <- clinical %>%
  left_join(dplyr::select(cluster_mapping, bio_id_merge, cluster), by = "bio_id_merge") %>%
  filter(cluster != "Intermediate") %>%
  dplyr::select(cluster, all_of(sig_diff_muts)) %>%
  pivot_longer(cols = -cluster, names_to = "gene", values_to = "mut") %>%
  filter(mut == 1) %>%
  summarize(mutation_event = n(), .by = c("gene", "cluster")) %>%
  left_join(cluster_sizes, by = "cluster") %>%
  mutate(frac = mutation_event / n * 100) %>%
  mutate(cluster = fct_relevel(cluster, c("Immature_like", "GMP_like", "Commited_like"))) %>%
  ggplot(aes(x = gene, y = frac, fill = cluster)) +
  geom_col(position = "dodge", col = "black") +
  scale_fill_manual(values = cluster_colors, labels = cluster_labels) +
  cowplot::theme_cowplot() +
  labs(x = "", y = "Fraction mutated (%)")

ggsave(file.path(fig_dir, "Fig2B_mutation_frequencies.pdf"), p_mut_freq,
       width = 6, height = 4, dpi = 300)

# =============================================================================
# 4. Mutations on diffusion map (Fig 2C)
# =============================================================================

p_mut_dm <- dm_coords %>%
  left_join(clinical, by = "bio_id_merge") %>%
  dplyr::select(DC1, DC2, FLT3_ITD_PCR, DNMT3A_R882, ras, PTPN11, RAD21) %>%
  pivot_longer(cols = FLT3_ITD_PCR:RAD21, names_to = "gene", values_to = "mut") %>%
  filter(!is.na(mut)) %>%
  mutate(mut = as.factor(mut)) %>%
  {
    ggplot(., aes(DC1, DC2)) +
      geom_point(data = filter(., mut == 0), col = "#dfe8dbff", alpha = 0.6) +
      geom_point(data = filter(., mut == 1), col = "#5cb350") +
      cowplot::theme_cowplot() +
      facet_wrap(. ~ gene)
  }

ggsave(file.path(fig_dir, "Fig2C_mutation_dm.pdf"), p_mut_dm,
       width = 7, height = 5, dpi = 300)

# =============================================================================
# 5. Differentiation scores by mutation (Fig 2D)
# =============================================================================

p_scores_by_mut <- diff_scores_wide %>%
  left_join(clinical, by = "bio_id_merge") %>%
  dplyr::select(Committed_like, Immature_like, FLT3_ITD_PCR, ras,
                DNMT3A_R882, PTPN11, RAD21) %>%
  pivot_longer(cols = FLT3_ITD_PCR:RAD21, names_to = "gene", values_to = "mutation_status") %>%
  pivot_longer(cols = c(Committed_like, Immature_like), names_to = "score_name",
               values_to = "score_value") %>%
  filter(!is.na(mutation_status)) %>%
  mutate(mutation_status = as.factor(mutation_status)) %>%
  ggplot(aes(mutation_status, score_value, fill = mutation_status)) +
  geom_boxplot() +
  facet_grid(score_name ~ gene, scales = "free") +
  ggpubr::stat_compare_means(comparisons = list(c("0", "1"))) +
  cowplot::theme_cowplot() +
  scale_fill_manual(values = c("#219ebc", "#e63946")) +
  theme(legend.position = "none") +
  labs(x = "", y = "Vector score")

ggsave(file.path(fig_dir, "Fig2D_scores_by_mutation.pdf"), p_scores_by_mut,
       width = 8, height = 5)

# =============================================================================
# 6. Survival by clinical features (Fig S2E-H)
# =============================================================================

surv_df <- cluster_mapping %>%
  dplyr::select(bio_id_merge, cluster) %>%
  left_join(clinical, by = "bio_id_merge")

# --- Fig S2E: BCL2 and HOX/Menin expression by cluster ---
bcl2_score <- as_tibble(vsn["BCL2", , drop = FALSE] %>% t(), rownames = "bio_id_merge") %>%
  dplyr::rename(BCL2_exp = BCL2)

menin_score <- as_tibble(
  t(vsn[grep("^HOXA|^HOXB|^MEIS1", rownames(vsn), value = TRUE), ]),
  rownames = "bio_id_merge"
) %>%
  pivot_longer(cols = -bio_id_merge) %>%
  summarize(Menin_score = mean(value), .by = bio_id_merge)

p_bcl2_menin <- cluster_mapping %>%
  left_join(bcl2_score) %>%
  left_join(menin_score) %>%
  mutate(cluster = factor(cluster, levels = c("Commited_like", "GMP_like", "Immature_like"))) %>%
  pivot_longer(cols = c(BCL2_exp, Menin_score), names_to = "score", values_to = "val") %>%
  filter(!is.na(cluster)) %>%
  ggplot(aes(cluster, val, fill = score)) +
  geom_boxplot() +
  scale_fill_manual(values = c("BCL2_exp" = "#5bb14fff", "Menin_score" = "#219ebc")) +
  cowplot::theme_cowplot() +
  labs(x = "", y = "Score")

ggsave(file.path(fig_dir, "FigS2E_bcl2_menin.pdf"), p_bcl2_menin,
       width = 6, height = 4, dpi = 300)

# --- Fig S2F: EFS by FLT3-ITD ---
pdf(file.path(fig_dir, "FigS2F_efs_flt3.pdf"), width = 5, height = 4.5)
ggsurvplot(survfit(Surv(efs_days, efsstat) ~ FLT3_ITD_PCR, surv_df),
           surv_df, pval = TRUE, risk.table = TRUE,
           palette = c("#219ebc", "#e63946"))
dev.off()

# --- Fig S2G: EFS by cohesin ---
pdf(file.path(fig_dir, "FigS2G_efs_cohesin.pdf"), width = 5, height = 4.5)
ggsurvplot(survfit(Surv(efs_days, efsstat) ~ cohesin, surv_df),
           surv_df, pval = TRUE, risk.table = TRUE,
           palette = c("#219ebc", "#e63946"))
dev.off()

# --- Fig S2H: EFS by ELN risk and triple-hit ---
pdf(file.path(fig_dir, "FigS2H_efs_eln.pdf"), width = 5, height = 4.5)
ggsurvplot(survfit(Surv(efs_days, efsstat) ~ ELN2022_risk, surv_df),
           surv_df, pval = TRUE, risk.table = TRUE,
           palette = c("#ffb703", "#219ebc", "#e63946"))
dev.off()

pdf(file.path(fig_dir, "FigS2H_efs_triple.pdf"), width = 5, height = 4.5)
surv_df %>%
  mutate(triple_status = case_when(
    FLT3_ITD_PCR == 1 & DNMT3A_R882 == 1 ~ "triple_hit",
    FLT3_ITD_PCR == 1 & DNMT3A_R882 == 0 ~ "FLT3_ITD",
    FLT3_ITD_PCR == 0 & DNMT3A_R882 == 1 ~ "DNMT3A",
    TRUE ~ "none"
  )) %>%
  filter(triple_status != "none") %>%
  ggsurvplot(survfit(Surv(efs_days, efsstat) ~ triple_status, .), .,
             pval = TRUE, risk.table = TRUE,
             palette = c("#ffb703", "#219ebc", "#e63946"))
dev.off()

cat("Done. Figures 2 and S2 saved.\n")
