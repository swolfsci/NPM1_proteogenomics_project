# =============================================================================
# 03_mito_score.R
#
# Mito-score computation, survival analysis, Cox models, discordant group
# characterization, and comparative GSEA.
# Generates: Figures 3 and S3
#
# Input:  output/preprocessed_data.RDS
#         output/differentiation_data.RDS
#         data/MitovsAll.csv
#         data/Human.MitoCarta3.0.xls
#         data/AMLCellType_Genesets.gmt
# Output: output/figures/Fig3_*.pdf
#         output/figures/FigS3_*.pdf
#         output/mito_score_data.RDS
# =============================================================================

library(tidyverse)
library(GSVA)
library(matrixStats)
library(survival)
library(survminer)
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(ComplexHeatmap)
library(forestmodel)
library(tidycmprsk)
library(ggsurvfit)
library(patchwork)

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
zeng_sets <- diff_data$zeng_sets

# =============================================================================
# 2. Compute mito-score
# =============================================================================

cat("Computing mito-score...\n")

# Load DEGs from Jayavelu, Wolf, Buettner et al. Cancer Cell 2022
mito_degs_raw <- read_csv(file.path(data_dir, "MitovsAll.csv"), show_col_types = FALSE)

# 1. DEGs: logFC > 0, adj.P.Val <= 0.05, present in our proteome
degs <- mito_degs_raw %>%
  filter(logFC > 0, adj.P.Val <= 0.05) %>%
  filter(PG.Genes %in% rownames(vsn)) %>%
  pull(PG.Genes)

# 2. Top quartile most variable proteins (HVPs)
protein_vars <- rowVars(vsn)
names(protein_vars) <- rownames(vsn)
n_q25 <- floor(nrow(vsn) * 0.25)
hvps <- names(sort(protein_vars, decreasing = TRUE))[1:n_q25]

# 3. Intersection = mito-score gene set
mito_genes <- intersect(degs, hvps)
cat(sprintf("Mito-score gene set: %d DEGs ∩ %d HVPs = %d genes\n",
            length(degs), n_q25, length(mito_genes)))

# 4. ssGSEA scoring
mito_ssgsea <- gsva(ssgseaParam(vsn, list(mito = mito_genes)), verbose = FALSE)
mito_score <- setNames(as.numeric(mito_ssgsea[1, ]), colnames(mito_ssgsea))

# Add to cluster mapping
cluster_mapping$mito_score_ssgsea <- mito_score[cluster_mapping$bio_id_merge]

# 5. Stratification at 75th percentile
q75 <- quantile(mito_score, 0.75)

# Merge clinical data for survival analyses (select only needed cols from
# cluster_mapping to avoid FLT3_ITD_PCR column collision with clinical)
surv_df <- cluster_mapping %>%
  dplyr::select(bio_id_merge, cluster, DC1, DC2, mito_score_ssgsea) %>%
  left_join(clinical, by = "bio_id_merge")

surv_df$high_mito <- surv_df$mito_score_ssgsea >= q75

# =============================================================================
# 3. Figure 3A: Mito-score on diffusion map
# =============================================================================

p_dm_mito <- cluster_mapping %>%
  ggplot(aes(DC1, DC2, col = mito_score_ssgsea)) +
  geom_point(size = 2) +
  cowplot::theme_cowplot() +
  scale_color_viridis_c(option = "inferno", name = "Mito-score")

ggsave(file.path(fig_dir, "Fig3A_dm_mito_score.pdf"), p_dm_mito,
       width = 6, height = 4)

# =============================================================================
# 4. Figure 3B: Mito-score by cluster
# =============================================================================

p_mito_cluster <- cluster_mapping %>%
  filter(cluster != "Intermediate") %>%
  ggplot(aes(x = cluster, y = mito_score_ssgsea, fill = cluster)) +
  geom_boxplot() +
  ggpubr::stat_compare_means(comparisons = list(
    c("Commited_like", "GMP_like"),
    c("Commited_like", "Immature_like"),
    c("GMP_like", "Immature_like")
  ), label = "p.signif") +
  scale_fill_manual(values = cluster_colors) +
  cowplot::theme_cowplot() +
  theme(legend.position = "none") +
  labs(x = "", y = "Mito-score (ssGSEA)")

ggsave(file.path(fig_dir, "Fig3B_mito_score_by_cluster.pdf"), p_mito_cluster,
       width = 5, height = 4)

# =============================================================================
# 5. Figures 3C-D: EFS by mito-score (all, FLT3-ITD negative)
# =============================================================================

p_efs_mito <- ggsurvplot(
  survfit(Surv(efs_months, efsstat) ~ high_mito, surv_df),
  surv_df, pval = TRUE, risk.table = TRUE,
  palette = c("#209ebc", "#E63A48")
)

p_efs_mito_flt3neg <- surv_df %>%
  filter(FLT3_ITD_PCR == 0) %>%
  {ggsurvplot(
    survfit(Surv(efs_months, efsstat) ~ high_mito, .),
    ., pval = TRUE, risk.table = TRUE,
    palette = c("#209ebc", "#E63A48")
  )}

pdf(file.path(fig_dir, "Fig3CD_efs_mito.pdf"), width = 10, height = 6)
arrange_ggsurvplots(list(p_efs_mito, p_efs_mito_flt3neg))
dev.off()

# =============================================================================
# 6. Figures S3A: OS by mito-score (all, FLT3-ITD negative)
# =============================================================================

p_os_mito <- ggsurvplot(
  survfit(Surv(os_months, stat) ~ high_mito, surv_df),
  surv_df, pval = TRUE, risk.table = TRUE,
  palette = c("#209ebc", "#E63A48")
)

p_os_mito_flt3neg <- surv_df %>%
  filter(FLT3_ITD_PCR == 0) %>%
  {ggsurvplot(
    survfit(Surv(os_months, stat) ~ high_mito, .),
    ., pval = TRUE, risk.table = TRUE,
    palette = c("#209ebc", "#E63A48")
  )}

pdf(file.path(fig_dir, "FigS3A_os_mito.pdf"), width = 10, height = 6)
arrange_ggsurvplots(list(p_os_mito, p_os_mito_flt3neg))
dev.off()

# =============================================================================
# 7. Figure S3B: Cox forest plots
# =============================================================================

# S3B: adjusted for FLT3 and cluster (binary)
pdf(file.path(fig_dir, "FigS3B_forest_efs_mito.pdf"), width = 6, height = 4)
surv_df %>%
  as.data.frame() %>%
  coxph(Surv(efs_days, efsstat) ~ high_mito + FLT3_ITD_PCR + cluster, .) %>%
  forest_model()
dev.off()

# S3C: extended model, FLT3-ITD negative, cluster adjustment
pdf(file.path(fig_dir, "FigS3C_forest_efs_mito_extended.pdf"), width = 6, height = 8)
surv_df %>%
  filter(FLT3_ITD_PCR == 0) %>%
  as.data.frame() %>%
  coxph(Surv(efs_days, efsstat) ~ high_mito + cluster + wbc_log10 +
          LDH_diag_log10 + age_cat + treatment_ITT + gender + type_AML, .) %>%
  forest_model()
dev.off()

# S3D: extended model with continuous differentiation scores
pdf(file.path(fig_dir, "FigS3D_forest_efs_mito_continuous.pdf"), width = 6, height = 8)
surv_df %>%
  left_join(diff_scores_wide, by = "bio_id_merge") %>%
  filter(FLT3_ITD_PCR == 0) %>%
  as.data.frame() %>%
  coxph(Surv(efs_days, efsstat) ~ high_mito + Immature_like + Committed_like +
          GMP_like + wbc_log10 + LDH_diag_log10 + age_cat + treatment_ITT +
          gender + type_AML, .) %>%
  forest_model()
dev.off()

# S3E: adjusted for FLT3 and continuous Immature_like
pdf(file.path(fig_dir, "FigS3E_forest_efs_mito_immature_cont.pdf"), width = 6, height = 4)
surv_df %>%
  left_join(diff_scores_wide, by = "bio_id_merge") %>%
  as.data.frame() %>%
  coxph(Surv(efs_days, efsstat) ~ high_mito + FLT3_ITD_PCR + Immature_like, .) %>%
  forest_model()
dev.off()

# =============================================================================
# 8. Figures S3F-G: Discordant group (non-immature, high mito-score)
# =============================================================================

surv_df <- surv_df %>%
  mutate(mito_group = case_when(
    high_mito & cluster == "Immature_like" ~ "high mito, immature",
    high_mito & cluster != "Immature_like" ~ "high mito, non-immature",
    TRUE ~ "rest"
  ))

# S3F: EFS by mito × immature groups
pdf(file.path(fig_dir, "FigS3F_efs_mito_immature_groups.pdf"), width = 6, height = 5)
ggsurvplot(
  survfit(Surv(efs_months, efsstat) ~ mito_group, surv_df),
  surv_df, pval = TRUE, risk.table = TRUE,
  palette = c("#023047", "#E63A48", "#5bb14fff")
)
dev.off()

# S3G: Discordant group on diffusion map
p_groups_dm <- surv_df %>%
  ggplot(aes(DC1, DC2, col = mito_group)) +
  geom_point(size = 2) +
  scale_color_manual(values = c(
    "high mito, immature" = "#023047",
    "high mito, non-immature" = "#E63A48",
    "rest" = "#5bb14fff"
  )) +
  cowplot::theme_cowplot() +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "FigS3G_dm_mito_groups.pdf"), p_groups_dm,
       width = 4.5, height = 4)

# =============================================================================
# 9. Figure S3H: Cumulative incidence of relapse
# =============================================================================

cir_df <- surv_df %>%
  mutate(
    cir_group = case_when(
      cluster != "Immature_like" & high_mito == TRUE  ~ "Non-immature, high Mito",
      cluster != "Immature_like" & high_mito == FALSE ~ "Non-immature, non-high Mito",
      cluster == "Immature_like" & high_mito == TRUE  ~ "Immature, high Mito",
      TRUE ~ NA_character_
    ),
    cuminc_factor = factor(cuminc, levels = c(0, 1, 2),
                           labels = c("Censored", "Relapse", "Death in CR"))
  )

cir_fit <- cuminc(Surv(efs_months, cuminc_factor) ~ cir_group, data = cir_df)

pdf(file.path(fig_dir, "FigS3H_cir_mito_groups.pdf"), width = 6, height = 4)
ggcuminc(cir_fit, outcome = "Relapse") +
  scale_color_manual(values = c(
    "Non-immature, high Mito"      = "#e63946",
    "Immature, high Mito"          = "#023047",
    "Non-immature, non-high Mito"  = "#5bb14f"
  )) +
  scale_fill_manual(values = c(
    "Non-immature, high Mito"      = "#e63946",
    "Immature, high Mito"          = "#023047",
    "Non-immature, non-high Mito"  = "#5bb14f"
  )) +
  labs(x = "Months", y = "Cumulative incidence of relapse") +
  add_risktable()
dev.off()

# =============================================================================
# 9b. Genetic context of the high-mito phenotype (Reviewer 2)
# =============================================================================
# Are high-mito cases enriched for any recurrent mutation, overall and within
# the FLT3-ITD-negative subset?

cat("Mito-score vs mutation enrichment...\n")

mutation_candidates <- intersect(
  c("FLT3_ITD_PCR", "FLT3_TKD", "ras", "DNMT3A_R882", "PTPN11", "RAD21",
    "SMC1A", "SMC3", "STAG2", "IDH1", "IDH2", "TET2", "WT1", "RUNX1",
    "ASXL1", "MYC", "NF1", "SRSF2", "cohesin"),
  colnames(clinical)
)

mito_mut <- cluster_mapping %>%
  mutate(high_mito = mito_score_ssgsea >= q75) %>%
  dplyr::select(bio_id_merge, high_mito, cluster) %>%
  inner_join(clinical, by = "bio_id_merge")

# genes with >= 3 mutated cases (Fisher on singletons is noise)
gene_cols_test <- mutation_candidates[
  map_lgl(mutation_candidates, ~ sum(mito_mut[[.x]] == 1, na.rm = TRUE) >= 3)
]

fisher_per_gene <- function(df) {
  map_dfr(gene_cols_test, function(g) {
    tab <- table(df[[g]], df$high_mito)
    if (any(dim(tab) < 2)) return(NULL)
    ft <- fisher.test(tab)
    tibble(gene = g,
           n_high   = sum(df[[g]] == 1 & df$high_mito, na.rm = TRUE),
           n_low    = sum(df[[g]] == 1 & !df$high_mito, na.rm = TRUE),
           frac_high = n_high / sum(df$high_mito, na.rm = TRUE),
           frac_low  = n_low  / sum(!df$high_mito, na.rm = TRUE),
           odds_ratio = unname(ft$estimate),
           p_value = ft$p.value)
  }) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    arrange(p_value)
}

if (length(gene_cols_test) > 0) {
  mito_mut_full <- fisher_per_gene(mito_mut)
  write_csv(mito_mut_full, file.path(table_dir, "mito_mutation_enrichment_full.csv"))

  if ("FLT3_ITD_PCR" %in% colnames(mito_mut)) {
    write_csv(fisher_per_gene(filter(mito_mut, FLT3_ITD_PCR == 0)),
              file.path(table_dir, "mito_mutation_enrichment_flt3neg.csv"))
  }

  p_mito_mut <- mito_mut_full %>%
    slice_min(p_value, n = 12) %>%
    pivot_longer(c(frac_high, frac_low), names_to = "group", values_to = "fraction") %>%
    mutate(group = recode(group, frac_high = "High Mito", frac_low = "Non-high Mito"),
           gene  = fct_reorder(gene, fraction, .fun = max, .desc = TRUE)) %>%
    ggplot(aes(gene, fraction, fill = group)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    scale_fill_manual(values = c("High Mito" = "#B2182B", "Non-high Mito" = "grey60")) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = NULL, y = "Fraction mutated", fill = NULL,
         title = "Mutation frequencies by Mito-AML status") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(fig_dir, "FigS3I_mito_mutation_enrichment.pdf"),
         p_mito_mut, width = 9, height = 5)
}

# =============================================================================
# 10. Figure 3E: Comparative GSEA — discordant group characterization
# =============================================================================

cat("Running comparative GSEA...\n")

# Define groups for GSEA
df_gsea <- cluster_mapping %>%
  dplyr::select(bio_id_merge, cluster, mito_score_ssgsea) %>%
  mutate(
    high_mito = mito_score_ssgsea >= q75,
    immature  = cluster == "Immature_like",
    group = case_when(
      !immature &  high_mito ~ "discordant",
      !immature & !high_mito ~ "reference",
       immature &  high_mito ~ "immature_high_mito",
       immature & !high_mito ~ "immature_non_high_mito"
    )
  ) %>%
  left_join(diff_scores_wide, by = "bio_id_merge") %>%
  left_join(clinical, by = "bio_id_merge")

# Prepare gene set databases
hallmark_df <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

mitocarta <- readxl::read_excel(file.path(data_dir, "Human.MitoCarta3.0.xls"),
                                sheet = 4) %>%
  janitor::clean_names()

mitocarta_df <- mitocarta %>%
  mutate(genes = str_split(genes, ",\\s*")) %>%
  unnest(genes) %>%
  dplyr::select(gs_name = mito_pathway, gene_symbol = genes)

zeng_df <- enframe(zeng_sets, name = "gs_name", value = "gene_symbol") %>%
  unnest(gene_symbol)

# --- Contrast A: discordant vs immature_high_mito ---
idx_a1 <- df_gsea$group %in% c("discordant", "immature_high_mito")
dep_a1 <- run_limma_contrast(vsn[, df_gsea$bio_id_merge[idx_a1]],
                              df_gsea[idx_a1, ], "group",
                              "discordant - immature_high_mito")
gsea_a1 <- run_gsea_all(dep_a1, hallmark_df, mitocarta_df, zeng_df)

# --- Contrast A2: discordant vs reference ---
idx_a2 <- df_gsea$group %in% c("discordant", "reference")
dep_a2 <- run_limma_contrast(vsn[, df_gsea$bio_id_merge[idx_a2]],
                              df_gsea[idx_a2, ], "group",
                              "discordant - reference")
gsea_a2 <- run_gsea_all(dep_a2, hallmark_df, mitocarta_df, zeng_df)

# --- Panel B: FLT3-ITD negative subgroup ---
df_gsea_neg <- df_gsea %>% filter(FLT3_ITD_PCR == 0)

idx_b1 <- df_gsea_neg$group %in% c("discordant", "immature_high_mito")
dep_b1 <- run_limma_contrast(vsn[, df_gsea_neg$bio_id_merge[idx_b1]],
                              df_gsea_neg[idx_b1, ], "group",
                              "discordant - immature_high_mito")
gsea_b1 <- run_gsea_all(dep_b1, hallmark_df, mitocarta_df, zeng_df)

idx_b2 <- df_gsea_neg$group %in% c("discordant", "reference")
dep_b2 <- run_limma_contrast(vsn[, df_gsea_neg$bio_id_merge[idx_b2]],
                              df_gsea_neg[idx_b2, ], "group",
                              "discordant - reference")
gsea_b2 <- run_gsea_all(dep_b2, hallmark_df, mitocarta_df, zeng_df)

# --- Panel C: Residualized mito-score ---
cat("Residualized mito-score GSEA...\n")
fit_resid <- lm(mito_score_ssgsea ~ Immature_like + Committed_like + GMP_like,
                data = df_gsea)
df_gsea$mito_resid <- residuals(fit_resid)
cat(sprintf("R² explained by differentiation: %.3f\n", summary(fit_resid)$r.squared))

design_c <- model.matrix(~ mito_resid, data = df_gsea)
fit_c <- lmFit(vsn[, df_gsea$bio_id_merge], design_c) %>% eBayes()
dep_c <- topTable(fit_c, coef = "mito_resid", number = Inf, sort.by = "none") %>%
  rownames_to_column("protein") %>% as_tibble()
gsea_c <- run_gsea_all(dep_c, hallmark_df, mitocarta_df, zeng_df)

# Select pathways
top_a <- select_paired_pathways(gsea_a1, gsea_a2)
top_b <- select_paired_pathways(gsea_b1, gsea_b2)
top_c <- select_top_pathways(gsea_c) %>% clean_labels()

cat(sprintf("Selected pathways: A=%d, B=%d, C=%d\n",
            nrow(top_a), nrow(top_b), nrow(top_c)))

# Group sizes for subtitles
n_disc_full <- sum(df_gsea$group == "discordant")
n_ref_full  <- sum(df_gsea$group == "reference")
n_imm_full  <- sum(df_gsea$group == "immature_high_mito")
n_disc_neg  <- sum(df_gsea$group == "discordant" & df_gsea$FLT3_ITD_PCR == 0)
n_ref_neg   <- sum(df_gsea$group == "reference" & df_gsea$FLT3_ITD_PCR == 0)
n_imm_neg   <- sum(df_gsea$group == "immature_high_mito" & df_gsea$FLT3_ITD_PCR == 0)

# --- Figure 3E: Paired dotplot, full cohort ---
p_gsea_a <- make_paired_dotplot(
  top_a, gsea_a1, gsea_a2,
  title = "High mito-score phenotype: comparative GSEA",
  subtitle = sprintf(
    "Non-immature high-mito (n=%d) vs immature high-mito (n=%d) and non-immature reference (n=%d)",
    n_disc_full, n_imm_full, n_ref_full
  )
)

ggsave(file.path(fig_dir, "Fig3E_mito_phenotype.pdf"), p_gsea_a,
       width = 6, height = 8)

# --- Figures S3K-L: FLT3-neg paired + residualized lollipop ---
p_gsea_b <- make_paired_dotplot(
  top_b, gsea_b1, gsea_b2,
  title = "FLT3-ITD negative subgroup",
  subtitle = sprintf(
    "Non-immature high-mito (n=%d) vs immature high-mito (n=%d) and reference (n=%d)",
    n_disc_neg, n_imm_neg, n_ref_neg
  )
)

p_gsea_c <- make_lollipop_plot(
  top_c,
  "Residualized mito-score",
  sprintf("Differentiation-independent component (n=%d)", nrow(df_gsea)),
  up_label = "Positive association",
  down_label = "Negative association"
)

p_suppl <- p_gsea_b / p_gsea_c +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave(file.path(fig_dir, "FigS3KL_mito_phenotype.pdf"), p_suppl,
       width = 7, height = 10)

# Save selected pathways
bind_rows(
  top_a %>% clean_labels() %>% mutate(panel = "A_full_cohort"),
  top_b %>% clean_labels() %>% mutate(panel = "B_flt3_neg"),
  top_c %>% mutate(panel = "C_residualized")
) %>%
  dplyr::select(panel, source, ID, label,
                any_of(c("NES_c1", "padj_c1", "NES_c2", "padj_c2",
                          "NES", "p.adjust"))) %>%
  write_csv(file.path(table_dir, "selected_pathways.csv"))

# =============================================================================
# 11. Save mito-score data for downstream scripts
# =============================================================================

saveRDS(list(
  mito_genes = mito_genes,
  mito_score = mito_score,
  q75 = q75,
  cluster_mapping = cluster_mapping,
  gsea_results = list(
    gsea_a1 = gsea_a1, gsea_a2 = gsea_a2,
    gsea_b1 = gsea_b1, gsea_b2 = gsea_b2,
    gsea_c = gsea_c
  ),
  hallmark_df = hallmark_df,
  mitocarta_df = mitocarta_df,
  zeng_df = zeng_df
), file.path(output_dir, "mito_score_data.RDS"))

cat("Done. Figures 3 and S3 saved.\n")
