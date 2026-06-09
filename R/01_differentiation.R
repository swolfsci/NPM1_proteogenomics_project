# =============================================================================
# 01_differentiation.R
#
# Proteomic differentiation landscape of NPM1-mutated AML.
# Generates: Figures 1 and S1
#
# Input:  output/preprocessed_data.RDS (from 00_preprocessing.R)
#         data/AMLCellType_Genesets.gmt
#         data/Furtwaengler_signatures.csv
#         data/cluster_mapping.csv
# Output: output/differentiation_data.RDS
#         output/figures/Fig1_*.pdf
#         output/figures/FigS1_*.pdf
# =============================================================================

library(tidyverse)
library(GSVA)
library(ComplexHeatmap)
library(destiny)
library(survival)
library(survminer)
library(PCGSE)
library(patchwork)

source("R/utils.R")

# =============================================================================
# 1. Load data
# =============================================================================

cat("Loading preprocessed data...\n")
preprocessed <- readRDS(file.path(output_dir, "preprocessed_data.RDS"))
vsn <- preprocessed$vsn_matrix
clinical <- preprocessed$clinical

cluster_mapping <- read.csv(file.path(data_dir, "cluster_mapping.csv"))

# =============================================================================
# 2. Gene set preparation
# =============================================================================

# Zeng AML cell type signatures
zeng_sets <- GSVA::readGMT(file.path(data_dir, "AMLCellType_Genesets.gmt"),
                            valueType = "list")
set_vector <- str_detect(names(zeng_sets), "Top100")
zeng_sets <- c(zeng_sets[1:5], zeng_sets[set_vector])

# Furtwaengler single-cell proteomics signatures (from DOI: 10.1126/science.adr8785)
furtwaengler_gmdp <- c("AZU1", "PRTN3", "ELANE", "CTSG", "MPO", "CALR")
furtwaengler_mdp <- c("LGALS1", "PLD4", "SAMHD1", "ANXA2", "COL11A1", "LYZ",
                       "CORO1A", "VIM", "HLA-DRB1")
furtwaengler_predc <- c("CORO1A", "LSP1")

zeng_sets$furtwaengler_gmdp <- furtwaengler_gmdp
zeng_sets$furtwaengler_mdp <- furtwaengler_mdp
zeng_sets$furtwaengler_predc <- furtwaengler_predc
zeng_sets$`GMP-like-Top100` <- NULL

# =============================================================================
# 3. GSVA enrichment scoring
# =============================================================================

cat("Running GSVA on cell type signatures...\n")
zeng_gsva <- gsva(gsvaParam(exprData = vsn, geneSets = zeng_sets, kcdf = "Gaussian"))

# Build long-format GSVA scores for selected signatures
zeng_gsva[
  c("Mono-like-Top100", "MLL_LSC_Somervaille2009_Down",
    "MLL_LSC_Somervaille2009_Up", "LSPC-Primed-Top100",
    "ProMono-like-Top100", "cDC-like-Top100",
    grep("furtwaengler", rownames(zeng_gsva), value = TRUE)),
] %>%
  as_tibble(rownames = "signature") %>%
  pivot_longer(cols = -signature, names_to = "bio_id_merge", values_to = "score_value") %>%
  pivot_wider(id_cols = bio_id_merge, names_from = signature, values_from = score_value) %>%
  column_to_rownames("bio_id_merge") %>%
  scale() %>%
  as_tibble(rownames = "bio_id_merge") %>%
  pivot_longer(cols = -bio_id_merge, names_to = "name", values_to = "value") -> zeng_gsva_long

# Aggregate into differentiation scores
diff_scores_wide <- aggregate_differentiation_scores(zeng_gsva_long)
diff_scores_long <- diff_scores_wide %>%
  pivot_longer(cols = -bio_id_merge, names_to = "differentiation", values_to = "score") %>%
  mutate(across(.cols = score, ~ scale(.)[, 1]), .by = differentiation)

# =============================================================================
# 4. PCA and PCGSE
# =============================================================================

cat("Running PCA and PCGSE...\n")
hvps <- matrixStats::rowVars(vsn) %>%
  sort(decreasing = TRUE) %>%
  head(2000) %>%
  names()

# PCA
pca_scores <- vsn[hvps, ] %>%
  t() %>%
  pcaMethods::pca(nPcs = 20) %>%
  pcaMethods::scores() %>%
  as_tibble(rownames = "bio_id_merge")

# Full proteome PCA for biplot
pca_full <- prcomp(t(vsn), center = TRUE, scale. = FALSE)
var_explained <- summary(pca_full)$importance[2, 1:5] * 100

# PCGSE
mtx <- table(stack(zeng_sets)[2:1])
binary_matrix <- as.matrix((mtx > 0) + 0)
common_genes <- intersect(hvps, colnames(binary_matrix))
data_input <- t(vsn[common_genes, ])
binary_matrix_aligned <- binary_matrix[, common_genes]
set_sizes <- rowSums(binary_matrix_aligned)
binary_matrix_clean <- binary_matrix_aligned[set_sizes > 1, ]
pcgse_result <- PCGSE::pcgse(data = data_input, gene.sets = binary_matrix_clean, pc.indexes = 1:2)

# --- Fig S1B: PCGSE heatmap ---
pcgse_stats <- pcgse_result$statistics %>%
  as_tibble(rownames = "signature") %>%
  dplyr::rename(PC1 = V1, PC2 = V2) %>%
  mutate(across(c(PC1, PC2), ~ . * -1))

pdf(file.path(fig_dir, "FigS1B_pcgse_heatmap.pdf"), width = 4.5, height = 4.5)
pcgse_stats %>%
  column_to_rownames("signature") %>%
  as.matrix() %>%
  Heatmap(col = circlize::colorRamp2(c(-3, 0, 3), c("#0b97baff", "#f1faee", "#d5510aff")))
dev.off()

# --- Fig S1C: PCGSE lollipop ---
pcgse_pvals <- pcgse_result$p.values %>%
  as_tibble(rownames = "signature") %>%
  dplyr::rename(PC1 = V1, PC2 = V2) %>%
  pivot_longer(cols = -signature, names_to = "PC", values_to = "pval")

pcgse_stats_long <- pcgse_stats %>%
  pivot_longer(cols = -signature, names_to = "PC", values_to = "statistics")

p_pcgse <- left_join(pcgse_pvals, pcgse_stats_long, by = c("signature", "PC")) %>%
  ggplot(aes(x = tidytext::reorder_within(signature, statistics, PC),
             xend = tidytext::reorder_within(signature, statistics, PC),
             yend = statistics)) +
  geom_segment(y = 0, linewidth = 1) +
  geom_point(aes(y = statistics, size = -log10(pval), col = -log10(pval))) +
  tidytext::scale_x_reordered() +
  cowplot::theme_cowplot() +
  coord_flip() +
  facet_wrap(~ PC, scale = "free") +
  labs(x = "", y = "")

ggsave(file.path(fig_dir, "FigS1C_pcgse_lollipop.pdf"), p_pcgse, width = 12, height = 4)

# =============================================================================
# 5. Diffusion map
# =============================================================================

cat("Computing diffusion map...\n")
dm <- vsn[hvps, ] %>%
  t() %>%
  pcaMethods::pca(nPcs = 20) %>%
  pcaMethods::scores() %>%
  DiffusionMap(k = 283)

dm_coords <- dm %>%
  as.data.frame() %>%
  as_tibble(rownames = "bio_id_merge")

# --- Fig S1A: GSVA heatmap ---
hm_color <- circlize::colorRamp2(c(-0.75, 0, 0.75), c("#d5510aff", "#f1faee", "#0b97baff"))

categ_color <- list(
  "FLT3_ITD_PCR" = c("1" = "#023047", "0" = "#fefae0"),
  "ELN2022_risk" = c("adverse" = "#e63946", "favorable" = "#219ebc", "intermediate" = "#fcbf49"),
  "gender" = c("female" = "#219ebc", "male" = "#ffb703"),
  "treatment_ITT" = c("Arm A: ATRA" = "#ae2012", "Arm B: GO + ATRA" = "#0a9396"),
  "normal" = c("1" = "#fefae0", "0" = "#023047"),
  "DNMT3A" = c("1" = "#023047", "0" = "#fefae0"),
  "type_AML" = c("deNovo AML" = "#fefae0", "sAML" = "#e63946", "tAML" = "#219ebc"),
  "age_cat" = c("<50" = "#ffb703", ">=50 & <70" = "#219ebc", ">=70" = "#023047")
)

hm_anno_df <- clinical %>%
  filter(bio_id_merge %in% colnames(zeng_gsva)) %>%
  arrange(match(bio_id_merge, colnames(zeng_gsva))) %>%
  dplyr::select(FLT3_ITD_PCR, ELN2022_risk, type_AML, age_cat, gender, normal,
                treatment_ITT, DNMT3A)

zeng_heatmap_top_anno <- HeatmapAnnotation(df = hm_anno_df, col = categ_color)

pdf(file.path(fig_dir, "FigS1A_gsva_heatmap.pdf"), width = 12, height = 6)
draw(Heatmap(zeng_gsva, col = hm_color, show_column_names = FALSE,
             top_annotation = zeng_heatmap_top_anno))
dev.off()

# --- Fig S1E: Diffusion maps colored by individual GSVA scores ---
p_dm_gsva <- dm_coords %>%
  left_join(zeng_gsva_long) %>%
  ggplot(aes(DC1, DC2, col = value)) +
  geom_point() +
  cowplot::theme_cowplot() +
  facet_wrap(. ~ name, nrow = 1) +
  scale_color_gradient2(high = "#d5510aff", mid = "#f1faee", low = "#0b97baff",
                        breaks = c(-3, 0, 2))

ggsave(file.path(fig_dir, "FigS1E_dm_gsva_signatures.pdf"), p_dm_gsva,
       width = 20, height = 3, dpi = 300)

# --- Fig 1C: Diffusion maps colored by aggregated differentiation scores ---
p_dm_diff <- dm_coords %>%
  left_join(diff_scores_long) %>%
  mutate(differentiation = fct_relevel(differentiation,
                                        c("Immature_like", "GMP_like", "Committed_like"))) %>%
  ggplot(aes(DC1, DC2, col = score)) +
  geom_point() +
  facet_wrap(. ~ differentiation) +
  cowplot::theme_cowplot() +
  scale_color_gradient2(high = "#d5510aff", mid = "#f1faee", low = "#0b97baff",
                        breaks = c(-3, 0, 2))

ggsave(file.path(fig_dir, "Fig1C_dm_differentiation.pdf"), p_dm_diff,
       width = 8, height = 3, dpi = 300)

# --- Fig 1D: Diffusion map with cluster assignments ---
p_dm_clusters <- dm_coords %>%
  left_join(dplyr::select(cluster_mapping, bio_id_merge, cluster), by = "bio_id_merge") %>%
  ggplot(aes(DC1, DC2, col = cluster)) +
  geom_point() +
  cowplot::theme_cowplot() +
  scale_color_manual(values = cluster_colors, labels = cluster_labels)

ggsave(file.path(fig_dir, "Fig1D_dm_clusters.pdf"), p_dm_clusters,
       width = 5, height = 3.5, dpi = 300)

# =============================================================================
# 6. Survival analysis by cluster
# =============================================================================

# --- Fig 1E: KM by cluster ---
surv_df <- cluster_mapping %>%
  dplyr::select(bio_id_merge, cluster) %>%
  left_join(clinical, by = "bio_id_merge")

pdf(file.path(fig_dir, "Fig1E_efs_by_cluster.pdf"), width = 6.5, height = 6)
ggsurvplot(
  survfit(Surv(efs_months, efsstat) ~ cluster, surv_df),
  surv_df,
  pval = TRUE, risk.table = TRUE,
  palette = unname(cluster_colors[c("Commited_like", "GMP_like", "Immature_like", "Intermediate")])
)
dev.off()

# --- Fig S2A: OS by cluster ---
pdf(file.path(fig_dir, "FigS2A_os_by_cluster.pdf"), width = 6.5, height = 6)
ggsurvplot(
  survfit(Surv(os_days, stat) ~ cluster, surv_df),
  surv_df,
  pval = TRUE, risk.table = TRUE,
  palette = unname(cluster_colors[c("Commited_like", "GMP_like", "Immature_like", "Intermediate")])
)
dev.off()

# --- Fig S2B: Multivariable Cox for EFS ---
cox_cluster <- surv_df %>%
  coxph(Surv(efs_days, efsstat) ~ cluster + ELN2022_risk + wbc_log10 +
          LDH_diag_log10 + age_cat + treatment_ITT + gender + type_AML, .)

p_forest_cluster <- forestmodel::forest_model(cox_cluster)
ggsave(file.path(fig_dir, "FigS2B_forest_cluster_efs.pdf"), p_forest_cluster,
       width = 6, height = 8)

# --- Fig S2C: Continuous differentiation scores Cox ---
cox_continuous <- diff_scores_wide %>%
  left_join(clinical) %>%
  coxph(Surv(efs_days, efsstat) ~ Committed_like + Immature_like + GMP_like, .)

p_forest_cont <- forestmodel::forest_model(cox_continuous)
ggsave(file.path(fig_dir, "FigS2C_forest_continuous.pdf"), p_forest_cont,
       width = 6, height = 4)

# --- Fig S2D: Immature vs rest ---
pdf(file.path(fig_dir, "FigS2D_efs_immature_vs_rest.pdf"), width = 5, height = 4.5)
surv_df %>%
  mutate(immature_like = cluster == "Immature_like") %>%
  ggsurvplot(
    survfit(Surv(efs_days, efsstat) ~ immature_like, .),
    ., pval = TRUE, risk.table = TRUE,
    palette = c("#219ebc", "#023047")
  )
dev.off()

# =============================================================================
# 7. PCA biplot (Fig S1H)
# =============================================================================

biplot_df <- tibble(
  bio_id_merge = rownames(pca_full$x),
  PC1 = pca_full$x[, 1],
  PC2 = pca_full$x[, 2]
) %>%
  left_join(dplyr::select(cluster_mapping, bio_id_merge, cluster), by = "bio_id_merge") %>%
  left_join(diff_scores_wide, by = "bio_id_merge")

score_names <- c("Immature_like", "Committed_like", "GMP_like")
cors <- map_dfr(score_names, function(s) {
  tibble(
    score = s,
    r_PC1 = cor(biplot_df$PC1, biplot_df[[s]], use = "complete.obs"),
    r_PC2 = cor(biplot_df$PC2, biplot_df[[s]], use = "complete.obs")
  )
})

arrow_scale <- max(abs(c(biplot_df$PC1, biplot_df$PC2))) * 0.7
arrows_df <- cors %>%
  mutate(
    x = r_PC1 * arrow_scale,
    y = r_PC2 * arrow_scale,
    label = str_replace(score, "_like", "-like") %>% str_replace("_", " ")
  )

p_biplot <- ggplot(biplot_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(col = cluster), size = 1.5, alpha = 0.7) +
  geom_segment(data = arrows_df,
               aes(x = 0, y = 0, xend = x, yend = y),
               arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
               linewidth = 0.9, color = "grey20") +
  geom_label(data = arrows_df,
             aes(x = x * 1.12, y = y * 1.12, label = label),
             size = 3.2, fontface = "bold", fill = "white",
             alpha = 0.85, linewidth = 0.2) +
  scale_color_manual(values = cluster_colors, labels = cluster_labels, name = "Cluster") +
  labs(
    x = sprintf("PC1 (%.1f%%)", var_explained[1]),
    y = sprintf("PC2 (%.1f%%)", var_explained[2]),
    title = "Proteome PCA with differentiation score projections"
  ) +
  coord_fixed() +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 12))

ggsave(file.path(fig_dir, "FigS1H_pca_biplot.pdf"), p_biplot, width = 5, height = 4)

# --- Fig S1G: GSVA scores stratified by cluster ---
p_gsva_by_cluster <- diff_scores_long %>%
  left_join(dplyr::select(cluster_mapping, bio_id_merge, cluster), by = "bio_id_merge") %>%
  ggplot(aes(cluster, score, fill = cluster)) +
  geom_boxplot() +
  facet_wrap(~ differentiation) +
  scale_fill_manual(values = cluster_colors) +
  cowplot::theme_cowplot() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(fig_dir, "FigS1G_gsva_by_cluster.pdf"), p_gsva_by_cluster,
       width = 10, height = 4)

# =============================================================================
# 8. Protein correlates of the differentiation axes (limma)
# =============================================================================
# Which proteins track each continuous phenotype? We regress the proteome on
# each differentiation score and on the two diffusion components individually,
# using the continuous score as the single covariate (this avoids the
# arbitrariness of the discrete clustering). Saved as a supplementary table.

cat("Protein-differentiation associations (limma)...\n")

fit_continuous_axis <- function(expr_tbl, coef_name) {
  tbl    <- expr_tbl %>% filter(bio_id_merge %in% colnames(vsn))
  expr   <- vsn[, tbl$bio_id_merge]
  design <- model.matrix(reformulate(coef_name), data = tbl)
  limma::lmFit(expr, design = design) %>%
    limma::eBayes() %>%
    limma::topTable(coef = coef_name, n = Inf) %>%
    as_tibble(rownames = "protein") %>%
    mutate(rank = sign(logFC) * -log10(adj.P.Val)) %>%
    arrange(desc(rank))
}

axis_inputs <- list(
  immature_like  = list(diff_scores_wide, "Immature_like"),
  committed_like = list(diff_scores_wide, "Committed_like"),
  gmp_like       = list(diff_scores_wide, "GMP_like"),
  dc1            = list(dplyr::select(cluster_mapping, bio_id_merge, DC1, DC2), "DC1"),
  dc2            = list(dplyr::select(cluster_mapping, bio_id_merge, DC1, DC2), "DC2")
)

tt_list <- imap(axis_inputs, ~ fit_continuous_axis(.x[[1]], .x[[2]]))

# The score-based and DC-based axes should agree (sanity check)
print(left_join(tt_list$committed_like, tt_list$dc1, by = "protein") %>%
        rstatix::cor_test(rank.x, rank.y, method = "spearman"))
print(left_join(tt_list$gmp_like, tt_list$dc2, by = "protein") %>%
        rstatix::cor_test(rank.x, rank.y, method = "spearman"))

writexl::write_xlsx(tt_list,
                    file.path(table_dir, "protein_differentiation_association.xlsx"))

# =============================================================================
# 9. Zeng et al. (Cancer Discov 2025) scHierarchy markers (optional)
# =============================================================================
# Independent differentiation-state marker sets projected onto the diffusion map
# as an external cross-validation of the map geometry. Requires the supplementary
# marker table; skipped with a message if it is not present.

sch_path <- file.path(data_dir, "zeng_scHierarchy_markers.xlsx")
if (file.exists(sch_path)) {
  cat("Zeng scHierarchy marker projection...\n")
  sch_markers <- as.list(readxl::read_excel(sch_path))
  sch_gsva <- gsva(gsvaParam(exprData = vsn, geneSets = sch_markers, kcdf = "Gaussian"))
  sch_long <- sch_gsva %>% as_tibble(rownames = "set") %>% pivot_longer(cols = -set)

  p_sch <- dm_coords %>%
    dplyr::select(bio_id_merge, DC1, DC2) %>%
    left_join(sch_long, by = c("bio_id_merge" = "name")) %>%
    ggplot(aes(DC1, DC2, col = value)) +
    geom_point() +
    facet_wrap(~ set) +
    scale_color_gradient2(high = "#d5510aff", mid = "#f1faee", low = "#0b97baff",
                          breaks = c(-3, 0, 2)) +
    cowplot::theme_cowplot()
  ggsave(file.path(fig_dir, "FigS1I_zeng_schierarchy_dm.pdf"), p_sch,
         width = 9, height = 6)
} else {
  cat("Zeng scHierarchy markers not found; skipping projection.\n")
}

# =============================================================================
# 10. Save intermediate objects
# =============================================================================

saveRDS(list(
  zeng_gsva = zeng_gsva,
  zeng_gsva_long = zeng_gsva_long,
  zeng_sets = zeng_sets,
  diff_scores_wide = diff_scores_wide,
  diff_scores_long = diff_scores_long,
  dm_coords = dm_coords,
  cluster_mapping = cluster_mapping,
  hvps = hvps
), file.path(output_dir, "differentiation_data.RDS"))

cat("Done. Saved: output/differentiation_data.RDS\n")
