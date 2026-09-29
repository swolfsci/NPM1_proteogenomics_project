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
library(rstatix)
library(forestmodel)

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
#
# Final (resubmission) version: rather than testing a hand-picked panel, we test
# every recurrent mutation (>= 5% of the cohort) against the three continuous
# differentiation scores with a Wilcoxon rank-sum test, BH-correct across the
# whole gene x score family, and only plot the genes that remain significant on
# at least one axis. This replaces the earlier uncorrected pairwise version.

# Candidate recurrent mutations: binary 0/1 columns present in the clinical table
mutation_candidates <- c("FLT3_ITD_PCR", "FLT3_TKD", "ras", "DNMT3A_R882",
                         "PTPN11", "RAD21", "SMC1A", "SMC3", "STAG2",
                         "IDH1", "IDH2", "TET2", "WT1", "RUNX1", "ASXL1",
                         "MYC", "NF1", "SRSF2", "cohesin")
mutation_candidates <- intersect(mutation_candidates, colnames(clinical))

# Frequency filter: keep mutations seen in >= 5% of the proteomic cohort
five_perc_mut <- cluster_mapping %>%
  dplyr::select(bio_id_merge) %>%
  left_join(clinical, by = "bio_id_merge") %>%
  dplyr::select(any_of(mutation_candidates)) %>%
  summarize(across(everything(), ~ mean(.x == 1, na.rm = TRUE) * 100)) %>%
  pivot_longer(everything(), names_to = "gene", values_to = "freq") %>%
  filter(freq >= 5) %>%
  pull(gene)

# Wilcoxon association of each mutation with the continuous scores (BH-adjusted)
geno_diff <- test_genotype_differentiation(diff_scores_wide, clinical, five_perc_mut)

# Genes significantly associated with at least one differentiation axis
sig_genes <- geno_diff %>%
  dplyr::select(gene, name, p.adj) %>%
  pivot_wider(names_from = name, values_from = p.adj) %>%
  filter(if_any(any_of(c("Committed_like", "Immature_like", "GMP_like")), ~ . < 0.05)) %>%
  pull(gene)

# Guard for dummy/example clinical data, where nothing may reach significance
if (length(sig_genes) == 0) {
  sig_genes <- head(intersect(c("FLT3_ITD_PCR", "DNMT3A_R882", "RAD21"), five_perc_mut), 3)
}

sig_vector <- geno_diff %>%
  filter(gene %in% sig_genes) %>%
  mutate(p.adj = round(p.adj, 3)) %>%
  rstatix::add_xy_position()

p_scores_by_mut <- diff_scores_wide %>%
  left_join(clinical, by = "bio_id_merge") %>%
  dplyr::select(Committed_like, Immature_like, GMP_like, all_of(sig_genes)) %>%
  pivot_longer(cols = all_of(sig_genes), names_to = "gene", values_to = "mutation_status") %>%
  pivot_longer(cols = c(Committed_like, Immature_like, GMP_like),
               names_to = "name", values_to = "value") %>%
  filter(!is.na(mutation_status)) %>%
  mutate(mutation_status = factor(ifelse(mutation_status == 1, "Mut", "WT"),
                                  levels = c("WT", "Mut"))) %>%
  ggplot(aes(mutation_status, value, fill = mutation_status)) +
  geom_boxplot() +
  facet_grid(name ~ gene, scales = "free") +
  cowplot::theme_cowplot() +
  scale_fill_manual(values = c("WT" = "#219ebc", "Mut" = "#e63946")) +
  theme(legend.position = "none") +
  labs(x = "", y = "Vector score") +
  ggpubr::stat_pvalue_manual(sig_vector, label = "p.adj.signif")

ggsave(file.path(fig_dir, "Fig2D_scores_by_mutation.pdf"), p_scores_by_mut,
       width = 10, height = 6)

# =============================================================================
# 5b. Genetic vs non-genetic determinants of differentiation (Reviewer 2)
# =============================================================================
# How much of the differentiation state (DC1/DC2) is captured by the recurrent
# mutations that individually associate with a score (`sig_genes`, above)?
# A bootstrap R^2 quantifies the genetic contribution; the residual is the
# non-genetic component.

if (length(sig_genes) >= 2 && requireNamespace("boot", quietly = TRUE)) {
  boot_data <- cluster_mapping %>%
    dplyr::select(bio_id_merge, DC1, DC2) %>%
    left_join(dplyr::select(clinical, bio_id_merge, all_of(sig_genes)),
              by = "bio_id_merge") %>%
    dplyr::select(-bio_id_merge) %>%
    drop_na()

  b_dc1 <- bootstrap_r2(dplyr::select(boot_data, -DC2), "DC1", R = 2000)
  b_dc2 <- bootstrap_r2(dplyr::select(boot_data, -DC1), "DC2", R = 2000)
  tibble(
    component = c("DC1", "DC2"),
    r2        = c(b_dc1$t0, b_dc2$t0),
    ci_low    = c(boot::boot.ci(b_dc1, type = "perc")$percent[4],
                  boot::boot.ci(b_dc2, type = "perc")$percent[4]),
    ci_high   = c(boot::boot.ci(b_dc1, type = "perc")$percent[5],
                  boot::boot.ci(b_dc2, type = "perc")$percent[5])
  ) %>%
    write_csv(file.path(table_dir, "variance_explained_by_mutations.csv"))
}

# --- Fig S2L: multivariable EFS forest for the GMP-like score ----------------
if (length(sig_genes) >= 1) {
  gmp_forest <- diff_scores_wide %>%
    left_join(clinical, by = "bio_id_merge") %>%
    dplyr::select(bio_id_merge, GMP_like, efs_days, efsstat, wbc_log10,
                  LDH_diag_log10, age_cat, treatment_ITT, gender, type_AML) %>%
    left_join(dplyr::select(clinical, bio_id_merge, all_of(sig_genes)),
              by = "bio_id_merge") %>%
    dplyr::select(-bio_id_merge) %>%
    mutate(across(all_of(sig_genes), as.factor)) %>%
    coxph(Surv(efs_days, efsstat) ~ ., .) %>%
    forest_model()
  ggsave(file.path(fig_dir, "FigS2L_gmp_forest.pdf"), gmp_forest, width = 5, height = 6)
}

# =============================================================================
# 5c. FAB morphology vs proteomic landscape (Fig 2E, Fig S2E-G)
# =============================================================================
# FAB class is compared with (i) the categorical clusters (Fig S2E) and (ii) the
# continuous scores (Fig S2F). The variance of DC1/DC2 is then partitioned
# between FAB morphology and the recurrent mutations with variancePartition
# (Fig 2E) and confirmed by nested linear models (Fig S2G). Skipped if the
# clinical table carries no FAB column.

fab_col <- intersect(c("FAB", "fab", "fab_klassifikation"), colnames(clinical))

if (length(fab_col) >= 1) {
  fab_col <- fab_col[1]

  FAB_LEVELS   <- c("M0/1", "M2", "M4", "M5", "M6")
  CLUST_LEVELS <- c("Immature_like", "GMP_like", "Intermediate", "Committed_like")

  fab_df <- cluster_mapping %>%
    dplyr::select(bio_id_merge, cluster, DC1, DC2) %>%
    left_join(diff_scores_wide, by = "bio_id_merge") %>%
    left_join(clinical, by = "bio_id_merge") %>%
    mutate(
      FAB = .data[[fab_col]],
      # numeric FAB codes -> "M2" etc.; pass through if already a string
      FAB = if (is.numeric(FAB)) paste0("M", FAB) else as.character(FAB),
      FAB = ifelse(FAB %in% c("M0", "M1"), "M0/1", FAB),
      cluster = ifelse(cluster == "Commited_like", "Committed_like", cluster)
    ) %>%
    mutate(
      FAB     = droplevels(factor(FAB, levels = FAB_LEVELS)),
      cluster = droplevels(factor(cluster, levels = CLUST_LEVELS))
    ) %>%
    filter(!is.na(FAB), !is.na(cluster))

  # --- 5c-i. FAB vs cluster (Fig S2E) ----------------------------------------
  # Global Fisher test, standardized residuals, and pairwise 2x2 Fisher tests
  # (one FAB class vs one cluster, BH-adjusted) for the odds ratios in the text.
  contingency_tab <- table(FAB = fab_df$FAB, Cluster = fab_df$cluster)
  set.seed(42)
  fisher_global <- fisher.test(contingency_tab, simulate.p.value = TRUE, B = 1e5)

  stdres_df <- as.data.frame.table(suppressWarnings(chisq.test(contingency_tab))$stdres) %>%
    rename(stdres = Freq) %>%
    mutate(obs = as.vector(contingency_tab),
           direction = case_when(stdres >  1.96 ~ "enriched",
                                 stdres < -1.96 ~ "depleted",
                                 TRUE           ~ "ns"))
  write_csv(stdres_df, file.path(table_dir, "fab_cluster_stdres.csv"))

  # both logical levels are forced so table() stays 2x2 when a cell is empty
  pairwise_fab_cluster <- expand_grid(FAB = levels(fab_df$FAB),
                                      Cluster = levels(fab_df$cluster)) %>%
    pmap_dfr(function(FAB, Cluster) {
      tab <- table(factor(fab_df$FAB == FAB, levels = c(FALSE, TRUE)),
                   factor(fab_df$cluster == Cluster, levels = c(FALSE, TRUE)))
      ft <- fisher.test(tab)
      tibble(FAB = FAB, Cluster = Cluster, n = tab["TRUE", "TRUE"],
             odds_ratio = unname(ft$estimate), ci_low = ft$conf.int[1],
             ci_high = ft$conf.int[2], p_value = ft$p.value)
    }) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH"))
  write_csv(pairwise_fab_cluster, file.path(table_dir, "fab_cluster_pairwise_fisher.csv"))

  prop_df <- as.data.frame.table(contingency_tab) %>%
    rename(FAB = FAB, Cluster = Cluster, n = Freq) %>%
    group_by(Cluster) %>%
    mutate(col_prop = n / sum(n)) %>%
    ungroup()
  write_csv(prop_df, file.path(table_dir, "fab_cluster_column_proportions.csv"))

  p_prop <- ggplot(prop_df, aes(Cluster, FAB, fill = col_prop)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.0f%%\n(n=%d)", 100 * col_prop, n)),
              size = 3.1, color = "grey15") +
    scale_fill_gradient(low = "#F7FBFF", high = "#08519C",
                        name = "Fraction\nof cluster",
                        labels = scales::percent_format(accuracy = 1),
                        limits = c(0, NA)) +
    labs(x = NULL, y = "FAB",
         title = "FAB composition of proteomic clusters",
         subtitle = sprintf("Column-normalized; Fisher's exact p = %.3g (simulated, B = 1e5)",
                            fisher_global$p.value)) +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(file.path(fig_dir, "FigS2E_fab_cluster_proportions.pdf"),
         p_prop, width = 5.5, height = 4.5)

  # --- 5c-ii. FAB vs continuous scores (Fig S2F) ------------------------------
  # Kruskal-Wallis (BH across the three scores) with epsilon^2; Dunn post-hoc
  # (BH within score) only for scores with a significant omnibus test.
  score_vars   <- c("Immature_like", "GMP_like", "Committed_like")
  score_labels <- c(Immature_like = "Immature-Score", GMP_like = "GMP-Score",
                    Committed_like = "Committed-Score")

  kw_results <- map_dfr(score_vars, function(sv) {
    res <- kruskal.test(reformulate("FAB", sv), data = fab_df)
    tibble(score = sv, statistic = unname(res$statistic),
           df = unname(res$parameter), p_value = res$p.value,
           epsilon_sq = kruskal_effsize(fab_df, reformulate("FAB", sv))$effsize)
  }) %>% mutate(p_adj = p.adjust(p_value, method = "BH"))
  write_csv(kw_results, file.path(table_dir, "fab_kruskal_wallis.csv"))

  sig_scores <- kw_results %>% filter(p_adj < 0.05) %>% pull(score)
  dunn_results <- if (length(sig_scores) > 0) {
    map_dfr(sig_scores, function(sv) {
      fab_df %>% dunn_test(reformulate("FAB", sv), p.adjust.method = "BH") %>%
        mutate(score = sv)
    })
  } else tibble()
  if (nrow(dunn_results) > 0)
    write_csv(dunn_results, file.path(table_dir, "fab_dunn_posthoc.csv"))

  fab_long <- fab_df %>%
    dplyr::select(FAB, all_of(score_vars)) %>%
    pivot_longer(all_of(score_vars), names_to = "score", values_to = "value") %>%
    mutate(score = factor(score, levels = score_vars, labels = score_labels))

  p_fab_scores <- ggplot(fab_long, aes(FAB, value, fill = FAB)) +
    geom_boxplot(outlier.size = 0.6, alpha = 0.85, linewidth = 0.3) +
    facet_wrap(~ score, nrow = 1, scales = "free_y") +
    scale_fill_brewer(palette = "Paired", guide = "none") +
    labs(x = "FAB subtype", y = "GSVA enrichment score",
         title = "Differentiation scores stratified by FAB",
         subtitle = sprintf("Kruskal-Wallis p_adj: %s",
                            paste(sprintf("%s = %.2g", score_labels[kw_results$score],
                                          kw_results$p_adj), collapse = "  |  "))) +
    theme_bw(base_size = 11) +
    theme(strip.background = element_rect(fill = "grey95", color = NA),
          panel.grid.minor = element_blank())

  # brackets for significant Dunn comparisons, stacked above each facet's max
  if (nrow(dunn_results) > 0) {
    dunn_sig <- dunn_results %>%
      filter(p.adj < 0.05) %>%
      mutate(score = factor(score, levels = score_vars, labels = score_labels)) %>%
      group_by(score) %>%
      mutate(y.position = max(fab_long$value[fab_long$score == score[1]]) +
               0.15 * row_number()) %>%
      ungroup()
    if (nrow(dunn_sig) > 0)
      p_fab_scores <- p_fab_scores +
        ggpubr::stat_pvalue_manual(dunn_sig, label = "p.adj.signif", tip.length = 0.01,
                                   size = 3, inherit.aes = FALSE)
  }
  ggsave(file.path(fig_dir, "FigS2F_fab_continuous_scores.pdf"),
         p_fab_scores, width = 8, height = 4.5)

  # --- 5c-iii. Morphology vs genetics: shared model data ----------------------
  gene_predictors <- intersect(
    c("DNMT3A_R882", "FLT3_ITD_PCR", "FLT3_TKD", "MYC", "NF1",
      "PTPN11", "RAD21", "SRSF2", "STAG2"),
    colnames(clinical)
  )

  df_complete <- NULL
  keep_gene   <- character(0)
  if (length(gene_predictors) >= 1) {
    df_complete <- fab_df %>%
      dplyr::select(FAB, DC1, DC2, all_of(gene_predictors)) %>%
      drop_na() %>%
      mutate(across(all_of(gene_predictors), ~ droplevels(factor(.x))),
             FAB = droplevels(FAB))

    # drop predictors whose minor class is too sparse to fit
    keep_gene <- gene_predictors[map_lgl(gene_predictors, function(g) {
      nlevels(df_complete[[g]]) >= 2 && min(table(df_complete[[g]])) >= 5
    })]
  }

  var_colors <- c("Morphology" = "#3B6E8F", "Genetics" = "#C5742E",
                  "Unattributed" = "grey75")

  variance_bar <- function(d, title, subtitle) {
    ggplot(d, aes(value, fct_relevel(component, c("DC2", "DC1")),
                  fill = factor(group, levels = names(var_colors)))) +
      geom_col(width = 0.55, color = "white", linewidth = 0.4) +
      geom_text(aes(label = sprintf("%.1f%%", 100 * value)),
                position = position_stack(vjust = 0.5),
                color = "white", fontface = "bold", size = 4) +
      scale_fill_manual(values = var_colors, name = NULL) +
      scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                         expand = expansion(mult = c(0, 0.02))) +
      labs(x = "Fraction of variance", y = NULL, title = title, subtitle = subtitle) +
      cowplot::theme_cowplot() +
      theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
            legend.position = "top")
  }

  # --- 5c-iv. variancePartition (Fig 2E) --------------------------------------
  # Linear mixed model per component with FAB and each mutation as random
  # effects; mutation terms are summed into a genetic component.
  if (length(keep_gene) > 0 && requireNamespace("variancePartition", quietly = TRUE)) {
    vp_predictors <- c("FAB", keep_gene)
    y_mat <- t(as.matrix(df_complete[, c("DC1", "DC2")]))
    vp_form <- as.formula(paste("~", paste0("(1|", vp_predictors, ")", collapse = " + ")))
    vp_fit <- tryCatch(
      variancePartition::fitExtractVarPartModel(y_mat, vp_form,
                                                dplyr::select(df_complete, all_of(vp_predictors))),
      error = function(e) { message("  variancePartition failed: ", e$message); NULL })

    if (!is.null(vp_fit)) {
      vp_df <- as.data.frame(vp_fit) %>%
        rownames_to_column("component") %>%
        pivot_longer(-component, names_to = "predictor", values_to = "value") %>%
        mutate(group = case_when(predictor == "FAB" ~ "Morphology",
                                 predictor == "Residuals" ~ "Unattributed",
                                 TRUE ~ "Genetics"))
      write_csv(vp_df, file.path(table_dir, "fab_genes_variancepartition_per_predictor.csv"))

      vp_grouped <- vp_df %>%
        group_by(component, group) %>%
        summarise(value = sum(value), .groups = "drop")
      write_csv(vp_grouped, file.path(table_dir, "fab_genes_variancepartition.csv"))

      p_vp <- variance_bar(vp_grouped,
                           "Variance of DC1/DC2: morphology vs genetics",
                           sprintf("variancePartition (LMM; FAB + %d mutations as random effects); n = %d",
                                   length(keep_gene), nrow(df_complete)))
      ggsave(file.path(fig_dir, "Fig2E_variancepartition.pdf"), p_vp, width = 6, height = 3.5)
    }
  } else if (length(keep_gene) > 0) {
    message("variancePartition not installed; skipping Fig 2E.")
  }

  # --- 5c-v. Nested R^2 confirmation (Fig S2G) --------------------------------
  # DC ~ FAB, DC ~ mutations, DC ~ FAB + mutations. Forward decomposition (FAB,
  # + mutations, unattributed) for the figure; the reverse increment (FAB given
  # mutations) and nested F-tests check each block's unique contribution.
  # 95% percentile CIs from 2,000 bootstrap resamples.
  if (length(keep_gene) > 0 && requireNamespace("boot", quietly = TRUE)) {
    nested_r2 <- function(d, dc) {
      r2 <- function(rhs) summary(lm(reformulate(rhs, dc), d))$r.squared
      fab <- r2("FAB"); genes <- r2(keep_gene); full <- r2(c("FAB", keep_gene))
      c(R2_FAB = fab, R2_genes = genes, R2_full = full,
        delta_genes_given_FAB = full - fab, delta_FAB_given_genes = full - genes,
        unattributed = 1 - full)
    }

    set.seed(42)
    r2_estimates <- map_dfr(c("DC1", "DC2"), function(dc) {
      bt <- boot::boot(df_complete, function(d, i) nested_r2(d[i, ], dc), R = 2000)
      map_dfr(seq_along(bt$t0), function(k) {
        ci <- boot::boot.ci(bt, type = "perc", index = k)$percent
        tibble(component = dc, quantity = names(bt$t0)[k], estimate = bt$t0[k],
               ci_low = ci[4], ci_high = ci[5])
      })
    })
    write_csv(r2_estimates, file.path(table_dir, "fab_genes_nested_r2.csv"))

    nested_ftests <- map_dfr(c("DC1", "DC2"), function(dc) {
      m <- function(rhs) lm(reformulate(rhs, dc), df_complete)
      m_fab <- m("FAB"); m_genes <- m(keep_gene); m_full <- m(c("FAB", keep_gene))
      tibble(component = dc,
             p_FAB_vs_null     = anova(m("1"), m_fab)[["Pr(>F)"]][2],
             p_genes_given_FAB = anova(m_fab, m_full)[["Pr(>F)"]][2],
             p_FAB_given_genes = anova(m_genes, m_full)[["Pr(>F)"]][2])
    })
    write_csv(nested_ftests, file.path(table_dir, "fab_genes_nested_r2_ftests.csv"))

    decomp_df <- r2_estimates %>%
      filter(quantity %in% c("R2_FAB", "delta_genes_given_FAB", "unattributed")) %>%
      transmute(component, value = estimate,
                group = dplyr::recode(quantity, R2_FAB = "Morphology",
                                      delta_genes_given_FAB = "Genetics",
                                      unattributed = "Unattributed"))

    p_decomp <- variance_bar(decomp_df,
                             "Nested R^2: morphology first, then mutations",
                             sprintf("lm(DC ~ FAB) vs lm(DC ~ FAB + %d mutations); n = %d",
                                     length(keep_gene), nrow(df_complete)))
    ggsave(file.path(fig_dir, "FigS2G_nested_r2.pdf"), p_decomp, width = 6, height = 3.5)
  }
}

# =============================================================================
# 5d. FLT3-ITD x differentiation interaction and discordant biology (Reviewer 2)
# =============================================================================

int_df <- cluster_mapping %>%
  dplyr::select(bio_id_merge, DC1, DC2) %>%
  left_join(clinical, by = "bio_id_merge") %>%
  mutate(across(c(DC1, DC2), ~ as.numeric(scale(.))))

cat("FLT3-ITD x (DC1 + DC2) interaction model:\n")
print(summary(coxph(Surv(efs_days, efsstat) ~ FLT3_ITD_PCR * (DC1 + DC2), int_df))$coefficients)

# Discordant biology: FLT3-ITD cases that nonetheless sit at the committed
# (high DC1) end of the map -> are these driven by a co-mutation?
itd_df <- cluster_mapping %>%
  dplyr::select(bio_id_merge, DC1, DC2) %>%
  left_join(clinical, by = "bio_id_merge") %>%
  filter(FLT3_ITD_PCR == 1)

if (nrow(itd_df) > 0) {
  high_dc1_ids <- itd_df %>%
    filter(DC1 >= quantile(DC1, 0.5, na.rm = TRUE)) %>%
    pull(bio_id_merge)

  como_genes <- intersect(mutation_candidates, colnames(itd_df))
  como_genes <- como_genes[map_lgl(como_genes, ~ sum(itd_df[[.x]], na.rm = TRUE) >= 5)]
  como_genes <- setdiff(como_genes, "FLT3_ITD_PCR")

  if (length(como_genes) > 0) {
    map_dfr(como_genes, function(g) {
      tab <- table(itd_df[[g]], itd_df$bio_id_merge %in% high_dc1_ids)
      if (any(dim(tab) < 2)) return(NULL)
      ft <- fisher.test(tab)
      tibble(gene = g, odds_ratio = unname(ft$estimate), p_value = ft$p.value)
    }) %>%
      mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
      arrange(p_value) %>%
      write_csv(file.path(table_dir, "flt3itd_committed_comutations.csv"))
  }
}

# =============================================================================
# 5e. FLT3-TKD as a distinct lesion (Reviewer, Fig S2)
# =============================================================================

if ("FLT3_TKD" %in% colnames(clinical)) {
  tkd_df <- cluster_mapping %>%
    dplyr::select(bio_id_merge, DC1, DC2) %>%
    left_join(diff_scores_wide, by = "bio_id_merge") %>%
    left_join(dplyr::select(clinical, bio_id_merge, FLT3_TKD), by = "bio_id_merge") %>%
    filter(!is.na(FLT3_TKD))

  # FLT3-TKD on the diffusion map
  p_tkd_dm <- tkd_df %>%
    {
      ggplot(., aes(DC1, DC2)) +
        geom_point(data = filter(., FLT3_TKD == 0), col = "#dfe8dbff", alpha = 0.6) +
        geom_point(data = filter(., FLT3_TKD == 1), col = "#5cb350") +
        cowplot::theme_cowplot()
    }
  ggsave(file.path(fig_dir, "FigS2_flt3tkd_dm.pdf"), p_tkd_dm, width = 4.5, height = 4)

  # FLT3-TKD vs differentiation scores
  p_tkd_scores <- tkd_df %>%
    pivot_longer(cols = c(Immature_like, GMP_like, Committed_like),
                 names_to = "score", values_to = "value") %>%
    mutate(FLT3_TKD = as.factor(FLT3_TKD)) %>%
    ggplot(aes(score, value, fill = FLT3_TKD)) +
    geom_boxplot() +
    cowplot::theme_cowplot() +
    scale_fill_manual(values = c("#219ebc", "#e63946")) +
    theme(legend.position = "none") +
    labs(x = "", y = "Vector score") +
    ggpubr::stat_compare_means()
  ggsave(file.path(fig_dir, "FigS2_flt3tkd_scores.pdf"), p_tkd_scores, width = 6, height = 4)

  # Prognostic impact of FLT3-TKD
  tkd_surv <- clinical %>% filter(!is.na(FLT3_TKD))
  cat("FLT3-TKD univariable EFS Cox:\n")
  print(summary(coxph(Surv(efs_days, efsstat) ~ FLT3_TKD, tkd_surv))$coefficients)
  if (all(c("os_days", "stat") %in% colnames(clinical))) {
    cat("FLT3-TKD univariable OS Cox:\n")
    print(summary(coxph(Surv(os_days, stat) ~ FLT3_TKD, tkd_surv))$coefficients)
  }
}

# =============================================================================
# 6. Survival by clinical features (Fig S2H-K)
# =============================================================================

surv_df <- cluster_mapping %>%
  dplyr::select(bio_id_merge, cluster) %>%
  left_join(clinical, by = "bio_id_merge")

# --- BCL2 and HOX/Menin expression by cluster (not in the current figure set) ---
bcl2_score <- as_tibble(vsn["BCL2", , drop = FALSE] %>% t(), rownames = "bio_id_merge") %>%
  dplyr::rename(BCL2_exp = BCL2)

menin_score <- as_tibble(
  t(vsn[grep("^HOXA|^HOXB|^MEIS1", rownames(vsn), value = TRUE), ]),
  rownames = "bio_id_merge"
) %>%
  pivot_longer(cols = -bio_id_merge) %>%
  summarize(Menin_score = mean(value), .by = bio_id_merge)

bcl2_menin_long <- cluster_mapping %>%
  left_join(bcl2_score) %>%
  left_join(menin_score) %>%
  mutate(cluster = factor(cluster, levels = c("Commited_like", "GMP_like", "Immature_like"))) %>%
  pivot_longer(cols = c(BCL2_exp, Menin_score), names_to = "score", values_to = "val") %>%
  filter(!is.na(cluster))

# Dunn post-hoc (BH-adjusted) per score, keep only significant pairs
dunn_results <- bcl2_menin_long %>%
  group_by(score) %>%
  dunn_test(val ~ cluster) %>%
  filter(p.adj.signif != "ns") %>%
  add_xy_position()

p_bcl2_menin <- bcl2_menin_long %>%
  ggplot(aes(cluster, val, fill = score)) +
  geom_boxplot() +
  scale_fill_manual(values = c("BCL2_exp" = "#5bb14fff", "Menin_score" = "#219ebc")) +
  cowplot::theme_cowplot() +
  labs(x = "", y = "Score") +
  ggpubr::stat_pvalue_manual(dunn_results)

ggsave(file.path(fig_dir, "FigS2_bcl2_menin.pdf"), p_bcl2_menin,
       width = 6, height = 4, dpi = 300)

# --- Fig S2H: EFS by FLT3-ITD ---
pdf(file.path(fig_dir, "FigS2H_efs_flt3.pdf"), width = 5, height = 4.5)
ggsurvplot(survfit(Surv(efs_days, efsstat) ~ FLT3_ITD_PCR, surv_df),
           surv_df, pval = TRUE, risk.table = TRUE,
           palette = c("#219ebc", "#e63946"))
dev.off()

# --- Fig S2I: EFS by cohesin ---
pdf(file.path(fig_dir, "FigS2I_efs_cohesin.pdf"), width = 5, height = 4.5)
ggsurvplot(survfit(Surv(efs_days, efsstat) ~ cohesin, surv_df),
           surv_df, pval = TRUE, risk.table = TRUE,
           palette = c("#219ebc", "#e63946"))
dev.off()

# --- Fig S2J-K: EFS by ELN risk and triple-hit ---
pdf(file.path(fig_dir, "FigS2J_efs_eln.pdf"), width = 5, height = 4.5)
ggsurvplot(survfit(Surv(efs_days, efsstat) ~ ELN2022_risk, surv_df),
           surv_df, pval = TRUE, risk.table = TRUE,
           palette = c("#ffb703", "#219ebc", "#e63946"))
dev.off()

pdf(file.path(fig_dir, "FigS2K_efs_triple.pdf"), width = 5, height = 4.5)
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
