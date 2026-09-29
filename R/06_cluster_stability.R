# =============================================================================
# 06_cluster_stability.R
#
# Stability of the locked differentiation metaclusters.
#
# The metaclusters were defined once (k-means, k = 8, on DC1/DC2, one random
# start, set.seed(123)), after which the eight k-means clusters were merged by
# hand into Immature-like, GMP-like, Committed-like and Intermediate. Every one of
# those choices is arbitrary. This script re-runs the clustering under
# alternative choices and asks two questions:
#
#   1. Membership: would the same patients end up in the same metacluster?
#   2. Outcome: does the GMP-like survival benefit depend on the exact partition?
#
# Perturbation families (each produces a patient partition):
#   seed       k-means k = 8, one start, 1,000 seeds (the locked recipe)
#   k          k-means k = 4..12 (without 8), one start, 100 seeds each
#   algorithm  k-means (100 starts), PAM, Ward, Gaussian mixture; k = 4..12
#   embedding  diffusion map rebuilt with 1,000-5,000 HVPs, 10-30 PCs,
#              global vs local kernel (k = 283 vs 50), 2 vs 3 DCs
#   subsample  500 x 80% patient subsamples, diffusion map rebuilt, then the
#              locked recipe (k-means k = 8, one start)
#
# Each partition is turned into metaclusters by two labelling rules:
#   anchored   each cluster takes the locked label held by most of its members
#              (tests whether the locked groups re-emerge; mimics labelling by eye)
#   score_tau  each cluster is named after the continuous differentiation score
#              with the highest cluster mean, if that mean is >= tau SD
#              (tau = 0.5, 0.75, 1); otherwise Intermediate. Uses no locked
#              labels, so it also tests the hand-labelling step.
#
# Input:  output/preprocessed_data.RDS (from 00_preprocessing.R)
#         output/differentiation_data.RDS (from 01_differentiation.R)
# Output: output/tables/cluster_stability_runs.csv
#         output/tables/cluster_stability_summary.csv
#         output/tables/cluster_stability_patient_consensus.csv
#         output/tables/cluster_stability_internal_validity.csv
#         output/tables/cluster_stability_consensus_outcome.csv
#         output/figures/FigS_cluster_stability.pdf
#         output/figures/FigS_cluster_stability_consensus_km.pdf
#
# Outcome columns need efs_days/efsstat, os_days/stat, age and ELN2022_risk in
# clinical_data.csv; without them only the membership analyses are run.
# =============================================================================

library(mclust)      # attached first so tidyverse masks its map()/count()
library(tidyverse)
library(survival)
library(survminer)
library(destiny)
library(cluster)
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
diff_scores_wide <- diff_data$diff_scores_wide

meta_levels <- c("Immature_like", "GMP_like", "Commited_like", "Intermediate")

cm <- cluster_mapping %>%
  mutate(cluster = recode(cluster, Committed_like = "Commited_like"))
ids <- cm$bio_id_merge
locked <- setNames(cm$cluster, ids)
X_locked <- as.matrix(dplyr::select(cm, DC1, DC2))
rownames(X_locked) <- ids
vsn <- vsn[, ids]

# Continuous scores, z-scored so the score_tau threshold is in SD units. Named
# with the metacluster spelling so argmax gives the label directly.
S <- diff_scores_wide %>%
  filter(bio_id_merge %in% ids) %>%
  arrange(match(bio_id_merge, ids)) %>%
  transmute(Immature_like = Immature_like, GMP_like = GMP_like,
            Commited_like = Committed_like) %>%
  mutate(across(everything(), ~ scale(.)[, 1])) %>%
  as.matrix()
rownames(S) <- ids

outcome_cols <- c("efs_days", "efsstat", "os_days", "stat", "age", "ELN2022_risk")
has_outcome <- all(outcome_cols %in% colnames(clinical))
if (!has_outcome) message("06_cluster_stability.R: outcome columns missing; ",
                          "running membership analyses only.")

surv_df <- tibble(bio_id_merge = ids) %>%
  left_join(dplyr::select(clinical, bio_id_merge, any_of(outcome_cols)),
            by = "bio_id_merge") %>%
  as.data.frame()
rownames(surv_df) <- ids

# =============================================================================
# 2. Helpers
# =============================================================================

# Diffusion map exactly as in 01_differentiation.R, with the knobs exposed
build_dm <- function(mat, n_hvp = 2000, n_pcs = 20, dm_k = ncol(mat) - 1, n_dc = 2) {
  hv <- names(sort(matrixStats::rowVars(mat, useNames = TRUE), decreasing = TRUE))[1:n_hvp]
  pcs <- pcaMethods::scores(pcaMethods::pca(t(mat[hv, ]), nPcs = n_pcs))
  dm <- suppressWarnings(DiffusionMap(pcs, k = min(dm_k, ncol(mat) - 1)))
  out <- as.matrix(as.data.frame(dm)[, paste0("DC", 1:n_dc)])
  rownames(out) <- colnames(mat)
  out
}

label_anchored <- function(part, ref) {
  maj <- tapply(ref, part, function(r) names(which.max(table(r))))
  unname(maj[as.character(part)])
}

label_score <- function(part, S, tau) {
  f <- factor(part)
  cen <- rowsum(S, f) / tabulate(f)
  lab <- ifelse(apply(cen, 1, max) >= tau, colnames(S)[max.col(cen, ties.method = "first")],
                "Intermediate")
  names(lab) <- levels(f)
  unname(lab[as.character(part)])
}

labellings <- c("anchored", "score_0.5", "score_0.75", "score_1")

apply_labelling <- function(rule, part, ids_run) {
  if (rule == "anchored") return(label_anchored(part, locked[ids_run]))
  label_score(part, S[ids_run, , drop = FALSE], as.numeric(sub("score_", "", rule)))
}

jaccard <- function(a, b) if (!any(a | b)) NA_real_ else sum(a & b) / sum(a | b)

# HR of one metacluster vs the rest
cox_hr <- function(formula, d) {
  fit <- tryCatch(coxph(formula, data = d), error = function(e) NULL)
  if (is.null(fit)) return(c(NA, NA, NA, NA))
  s <- summary(fit)
  c(s$conf.int["grp", c(1, 3, 4)], s$coefficients["grp", "Pr(>|z|)"])
}

outcome_row <- function(lab, ids_run) {
  empty <- tibble(efs_hr = NA, efs_lo = NA, efs_hi = NA, efs_p = NA,
                  efs_adj_hr = NA, efs_adj_p = NA, os_hr = NA, os_p = NA,
                  imm_efs_hr = NA, com_efs_hr = NA)
  if (!has_outcome) return(empty)
  d <- surv_df[ids_run, ]
  hr_for <- function(level, f) {
    d$grp <- as.numeric(lab == level)
    if (sum(d$grp) < 10 || sum(d$grp) > nrow(d) - 10) return(c(NA, NA, NA, NA))
    cox_hr(f, d)
  }
  e  <- hr_for("GMP_like", Surv(efs_days, efsstat) ~ grp)
  ea <- hr_for("GMP_like", Surv(efs_days, efsstat) ~ grp + age + ELN2022_risk)
  o  <- hr_for("GMP_like", Surv(os_days, stat) ~ grp)
  tibble(efs_hr = e[1], efs_lo = e[2], efs_hi = e[3], efs_p = e[4],
         efs_adj_hr = ea[1], efs_adj_p = ea[4], os_hr = o[1], os_p = o[4],
         imm_efs_hr = hr_for("Immature_like", Surv(efs_days, efsstat) ~ grp)[1],
         com_efs_hr = hr_for("Commited_like", Surv(efs_days, efsstat) ~ grp)[1])
}

# Summarise one partition under every labelling rule. GMP-like membership of
# full-cohort runs is kept for the per-patient consensus.
gmp_membership <- list()

evaluate_partition <- function(part, family, params, ids_run = ids) {
  map_dfr(labellings, function(rule) {
    lab <- apply_labelling(rule, part, ids_run)
    ref <- locked[ids_run]
    if (length(ids_run) == length(ids)) {
      gmp_membership[[length(gmp_membership) + 1]] <<-
        list(family = family, labelling = rule, gmp = lab == "GMP_like")
    }
    jac <- map_dbl(meta_levels, ~ jaccard(lab == .x, ref == .x))
    tibble(family = family, params = params, labelling = rule,
           n_patients = length(ids_run), n_clusters = n_distinct(part),
           ari_vs_locked = mclust::adjustedRandIndex(lab, ref),
           identical_to_locked = all(lab == ref),
           n_immature = sum(lab == "Immature_like"), n_gmp = sum(lab == "GMP_like"),
           n_committed = sum(lab == "Commited_like"),
           n_intermediate = sum(lab == "Intermediate"),
           jac_immature = jac[1], jac_gmp = jac[2], jac_committed = jac[3],
           jac_intermediate = jac[4]) %>%
      bind_cols(outcome_row(lab, ids_run))
  })
}

# =============================================================================
# 3. Locked reference and internal validity of the DC1/DC2 space
# =============================================================================

cat("Locked reference...\n")
locked_row <- evaluate_partition(match(locked, meta_levels), "locked", "as published")
gmp_membership <- list()  # the locked partition does not count towards consensus

# Does DC1/DC2 carry discrete cluster structure? Average silhouette and gap
# statistic for k-means, plus the number of components preferred by BIC.
cat("Internal validity (silhouette, gap, BIC)...\n")
d_locked <- dist(X_locked)
set.seed(1)
gap_stat <- clusGap(X_locked, FUNcluster = function(x, k) kmeans(x, k, nstart = 25, iter.max = 50),
               K.max = 12, B = 200, verbose = FALSE)
internal <- tibble(
  k = 2:12,
  silhouette = map_dbl(2:12, function(k) {
    set.seed(1)
    mean(silhouette(kmeans(X_locked, k, nstart = 50)$cluster, d_locked)[, "sil_width"])
  }),
  gap = gap_stat$Tab[2:12, "gap"], gap_se = gap_stat$Tab[2:12, "SE.sim"]
)
gap_k <- maxSE(gap_stat$Tab[, "gap"], gap_stat$Tab[, "SE.sim"], method = "firstSEmax")
bic <- mclust::Mclust(X_locked, G = 1:12, verbose = FALSE)
internal <- internal %>%
  mutate(gap_selected_k = gap_k, bic_selected_G = bic$G, bic_model = bic$modelName)
write_csv(internal, file.path(table_dir, "cluster_stability_internal_validity.csv"))
cat(sprintf("  best silhouette k = %d; gap selects k = %d; BIC selects G = %d (%s)\n",
            internal$k[which.max(internal$silhouette)], gap_k, bic$G, bic$modelName))

# =============================================================================
# 4. Perturbation families
# =============================================================================

runs <- list()

cat("Family 1: k-means seeds (k = 8, one start)...\n")
runs$seed <- map_dfr(1:1000, function(s) {
  set.seed(s)
  evaluate_partition(kmeans(X_locked, 8)$cluster, "seed", sprintf("k=8 seed=%d", s))
})

cat("Family 2: number of clusters (k = 4..12, one start)...\n")
runs$k <- map_dfr(setdiff(4:12, 8), function(k) map_dfr(1:100, function(s) {
  set.seed(s)
  evaluate_partition(kmeans(X_locked, k)$cluster, "k", sprintf("k=%d seed=%d", k, s))
}))

cat("Family 3: clustering algorithm (k = 4..12)...\n")
hc_ward <- hclust(d_locked, method = "ward.D2")
runs$algorithm <- map_dfr(4:12, function(k) {
  set.seed(1)
  parts <- list(
    kmeans_best = kmeans(X_locked, k, nstart = 100)$cluster,
    pam         = pam(X_locked, k, cluster.only = TRUE),
    ward        = cutree(hc_ward, k),
    gmm         = mclust::Mclust(X_locked, G = k, modelNames = "VVV", verbose = FALSE)$classification
  )
  imap_dfr(parts, ~ evaluate_partition(.x, "algorithm", sprintf("%s k=%d", .y, k)))
})

cat("Family 4: diffusion map construction...\n")
emb_grid <- expand_grid(n_hvp = c(1000, 2000, 3000, 5000), n_pcs = c(10, 20, 30),
                        dm_k = c(length(ids) - 1, 50), n_dc = 2) %>%
  bind_rows(tibble(n_hvp = 2000, n_pcs = 20, dm_k = length(ids) - 1, n_dc = 3))
runs$embedding <- pmap_dfr(emb_grid, function(n_hvp, n_pcs, dm_k, n_dc) {
  emb <- build_dm(vsn, n_hvp, n_pcs, dm_k, n_dc)
  set.seed(1)
  evaluate_partition(kmeans(emb, 8, nstart = 100)$cluster, "embedding",
                     sprintf("hvp=%d pcs=%d dm_k=%d dcs=%d", n_hvp, n_pcs, dm_k, n_dc))
})

cat("Family 5: 80% patient subsamples (diffusion map rebuilt)...\n")
runs$subsample <- map_dfr(1:500, function(b) {
  set.seed(b)
  sub <- sort(sample(ids, round(0.8 * length(ids))))
  emb <- build_dm(vsn[, sub])
  evaluate_partition(kmeans(emb, 8)$cluster, "subsample", sprintf("b=%d", b), ids_run = sub)
})

family_levels <- c("locked", "seed", "k", "algorithm", "embedding", "subsample")
all_runs <- bind_rows(locked_row, bind_rows(runs)) %>%
  mutate(family = factor(family, levels = family_levels))
write_csv(all_runs, file.path(table_dir, "cluster_stability_runs.csv"))

# =============================================================================
# 5. Summaries
# =============================================================================

locked_hr <- locked_row$efs_hr[1]

summary_tbl <- all_runs %>%
  filter(family != "locked") %>%
  group_by(family, labelling) %>%
  summarise(
    n_runs = n(),
    pct_identical_to_locked = 100 * mean(identical_to_locked),
    ari_median = median(ari_vs_locked),
    jac_gmp_median = median(jac_gmp, na.rm = TRUE),
    jac_gmp_q10 = quantile(jac_gmp, 0.1, na.rm = TRUE),
    jac_immature_median = median(jac_immature, na.rm = TRUE),
    jac_committed_median = median(jac_committed, na.rm = TRUE),
    jac_intermediate_median = median(jac_intermediate, na.rm = TRUE),
    n_gmp_median = median(n_gmp),
    pct_no_gmp_group = 100 * mean(n_gmp == 0),
    efs_hr_median = median(efs_hr, na.rm = TRUE),
    efs_hr_q10 = quantile(efs_hr, 0.1, na.rm = TRUE),
    efs_hr_q90 = quantile(efs_hr, 0.9, na.rm = TRUE),
    pct_efs_hr_below_1 = 100 * mean(efs_hr < 1, na.rm = TRUE),
    pct_efs_p_below_05 = 100 * mean(efs_hr < 1 & efs_p < 0.05, na.rm = TRUE),
    pct_efs_adj_p_below_05 = 100 * mean(efs_adj_hr < 1 & efs_adj_p < 0.05, na.rm = TRUE),
    # share of runs with a GMP-like HR further from 1 than the locked one
    pct_more_extreme_than_locked = 100 * mean(efs_hr <= locked_hr, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(summary_tbl, file.path(table_dir, "cluster_stability_summary.csv"))
print(as.data.frame(summary_tbl %>% filter(labelling %in% c("anchored", "score_0.75")) %>%
  dplyr::select(family, labelling, pct_identical_to_locked, ari_median, jac_gmp_median,
                n_gmp_median, efs_hr_median, pct_efs_p_below_05,
                pct_more_extreme_than_locked)), digits = 3)

# =============================================================================
# 6. Per-patient consensus
# =============================================================================
# Share of full-cohort runs in which a patient is GMP-like, averaged over the
# four full-cohort families so the 1,000 seed runs do not dominate.

cat("Per-patient consensus...\n")
consensus <- map_dfr(labellings, function(rule) {
  keep <- keep(gmp_membership, ~ .x$labelling == rule)
  fam <- map_chr(keep, "family")
  mem <- do.call(cbind, map(keep, "gmp"))
  per_family <- sapply(unique(fam), function(f) rowMeans(mem[, fam == f, drop = FALSE]))
  tibble(bio_id_merge = ids, labelling = rule, locked_cluster = unname(locked),
         gmp_frequency = rowMeans(per_family))
}) %>%
  mutate(consensus_group = cut(gmp_frequency, c(-Inf, 0.2, 0.8, Inf),
                               labels = c("rarely GMP-like (<20%)", "boundary (20-80%)",
                                          "core GMP-like (>=80%)")))
write_csv(consensus, file.path(table_dir, "cluster_stability_patient_consensus.csv"))
print(consensus %>% count(labelling, locked_gmp = locked_cluster == "GMP_like", consensus_group) %>%
        pivot_wider(names_from = consensus_group, values_from = n, values_fill = 0) %>%
        as.data.frame())

cons_anch <- consensus %>% filter(labelling == "anchored")

if (has_outcome) {
  cons_df <- surv_df %>% left_join(cons_anch, by = "bio_id_merge") %>%
    mutate(efs_months = efs_days / 30.44, grp = consensus_group)
  tidy_cox <- function(fit, model) {
    s <- summary(fit)
    tibble(model = model, term = rownames(s$conf.int), hr = s$conf.int[, 1],
           lo = s$conf.int[, 3], hi = s$conf.int[, 4], p = s$coefficients[, "Pr(>|z|)"])
  }
  consensus_outcome <- bind_rows(
    tidy_cox(coxph(Surv(efs_days, efsstat) ~ relevel(grp, "rarely GMP-like (<20%)"), cons_df),
             "EFS ~ consensus group"),
    tidy_cox(coxph(Surv(efs_days, efsstat) ~ I(gmp_frequency * 10), cons_df),
             "EFS ~ GMP-like frequency (per 10%)"),
    tidy_cox(coxph(Surv(efs_days, efsstat) ~ I(gmp_frequency * 10) + age + ELN2022_risk,
                   cons_df), "EFS ~ GMP-like frequency (per 10%) + age + ELN2022"),
    tidy_cox(coxph(Surv(efs_days, efsstat) ~ I(locked_cluster == "GMP_like"),
                   filter(cons_df, consensus_group != "boundary (20-80%)")),
             "EFS ~ locked GMP-like, boundary patients removed")
  )
  write_csv(consensus_outcome, file.path(table_dir, "cluster_stability_consensus_outcome.csv"))
  print(as.data.frame(consensus_outcome), digits = 3)

  pdf(file.path(fig_dir, "FigS_cluster_stability_consensus_km.pdf"), width = 6.5, height = 6)
  print(ggsurvplot(survfit(Surv(efs_months, efsstat) ~ consensus_group, cons_df),
                   data = cons_df, risk.table = TRUE, pval = TRUE, xlab = "Months",
                   ylab = "Event-free survival", legend.title = "",
                   legend.labs = levels(cons_df$consensus_group),
                   palette = c("grey60", "#fcd26b", "#ffb703")))
  dev.off()
}

# =============================================================================
# 7. Figure
# =============================================================================

family_labels <- c(seed = "Seed\n(k = 8)", k = "k = 4-12", algorithm = "Algorithm",
                   embedding = "Embedding", subsample = "80% subsample")
labelling_labels <- c(anchored = "Anchored to locked labels", `score_0.5` = "Score rule, tau = 0.5",
                      `score_0.75` = "Score rule, tau = 0.75", score_1 = "Score rule, tau = 1")

# A. Where the GMP-like calls are stable
p_a <- cm %>%
  left_join(cons_anch, by = "bio_id_merge") %>%
  arrange(gmp_frequency) %>%
  ggplot(aes(DC1, DC2)) +
  geom_point(aes(fill = gmp_frequency), shape = 21, size = 2.2, stroke = 0.2, colour = "white") +
  geom_point(data = ~ filter(.x, cluster == "GMP_like"), shape = 21, size = 2.9,
             stroke = 0.5, colour = "black", fill = NA) +
  scale_fill_gradient(low = "#e8e8e8", high = "#b37d00", limits = c(0, 1),
                      labels = scales::percent, name = "Runs called\nGMP-like") +
  labs(subtitle = "Black ring: locked GMP-like") +
  cowplot::theme_cowplot(10)

# B. Membership agreement with the locked metaclusters (anchored labelling)
p_b <- all_runs %>%
  filter(family != "locked", labelling == "anchored") %>%
  dplyr::select(family, starts_with("jac_")) %>%
  pivot_longer(-family, names_to = "metacluster", values_to = "jaccard") %>%
  mutate(metacluster = factor(cluster_labels[c(jac_immature = "Immature_like",
                                               jac_gmp = "GMP_like",
                                               jac_committed = "Commited_like",
                                               jac_intermediate = "Intermediate")[metacluster]],
                              levels = cluster_labels[meta_levels])) %>%
  ggplot(aes(family, jaccard, fill = metacluster)) +
  geom_boxplot(outlier.size = 0.4, linewidth = 0.3, position = position_dodge2(padding = 0.2)) +
  scale_fill_manual(values = setNames(cluster_colors[meta_levels], cluster_labels[meta_levels]),
                    name = NULL) +
  scale_x_discrete(labels = family_labels) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = "Jaccard overlap with locked group") +
  cowplot::theme_cowplot(10) + theme(legend.position = "top")

# C. Silhouette: no preferred k
p_c <- internal %>%
  ggplot(aes(k, silhouette)) +
  geom_line(linewidth = 0.5) + geom_point(size = 2) +
  geom_vline(xintercept = 8, linetype = 2, colour = "grey50") +
  scale_x_continuous(breaks = 2:12) +
  labs(y = "Mean silhouette width", x = "k (k-means on DC1/DC2)",
       subtitle = sprintf("Gap statistic selects k = %d; mixture BIC selects G = %d",
                          gap_k, bic$G)) +
  cowplot::theme_cowplot(10)

p_top <- (p_a | p_c) + plot_layout(widths = c(1.2, 1))

if (has_outcome) {
  # D. GMP-like EFS hazard ratio across all runs
  p_d <- all_runs %>%
    filter(family != "locked", !is.na(efs_hr)) %>%
    mutate(labelling = factor(labelling_labels[labelling], levels = labelling_labels)) %>%
    ggplot(aes(family, efs_hr)) +
    geom_hline(yintercept = 1, colour = "grey40", linewidth = 0.4) +
    geom_hline(yintercept = locked_hr, linetype = 2, colour = "#b37d00", linewidth = 0.5) +
    geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.3, fill = "#ffe3a3") +
    geom_jitter(data = ~ filter(.x, family %in% c("algorithm", "embedding")),
                width = 0.12, size = 0.8, alpha = 0.7) +
    facet_wrap(~ labelling, nrow = 1) +
    scale_y_log10() +
    scale_x_discrete(labels = family_labels) +
    labs(x = NULL, y = "EFS hazard ratio, GMP-like vs rest (log scale)",
         subtitle = sprintf("Dashed line: locked partition (HR %.2f)", locked_hr)) +
    cowplot::theme_cowplot(10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background = element_blank())
  p_all <- p_top / p_b / p_d + plot_annotation(tag_levels = "A") +
    plot_layout(heights = c(1, 0.9, 1))
  fig_h <- 13
} else {
  p_all <- p_top / p_b + plot_annotation(tag_levels = "A")
  fig_h <- 8.5
}

ggsave(file.path(fig_dir, "FigS_cluster_stability.pdf"), p_all, width = 11, height = fig_h)
ggsave(file.path(fig_dir, "FigS_cluster_stability.png"), p_all, width = 11, height = fig_h,
       dpi = 200)

cat("Done.\n")
