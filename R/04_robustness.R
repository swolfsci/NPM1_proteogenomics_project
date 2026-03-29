# =============================================================================
# 04_robustness.R
#
# Mito-score robustness analyses:
#   A. Cutoff sensitivity (50–85% quantile)
#   B. Gene set permutation test (2,000 random gene sets)
#   C. Bootstrap confidence intervals (1,000 patient resamples)
#   D. DEG filter sensitivity (adj.P threshold × logFC threshold grid)
#   E. HVP threshold / scoring method comparison
#   F. Continuous Cox regression (HR per SD, no cutoff)
#   G. Leave-one-out protein stability
#   H. Split-half cross-validation (500 random splits)
#
# Input:  output/preprocessed_data.RDS
#         output/differentiation_data.RDS
#         output/mito_score_data.RDS
#         data/MitovsAll.csv
# Output: output/tables/robustness_*.csv
#         output/robustness_results.RDS
#         output/figures/FigS3B_permutation.pdf
#         output/figures/FigS3C_cutoff_sensitivity.pdf
#         output/figures/TableS2.xlsx
# =============================================================================

library(tidyverse)
library(survival)
library(matrixStats)
library(GSVA)

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

mito_data <- readRDS(file.path(output_dir, "mito_score_data.RDS"))
mito_genes <- mito_data$mito_genes
obs_score <- mito_data$mito_score

mitovsall <- read_csv(file.path(data_dir, "MitovsAll.csv"), show_col_types = FALSE)

# Base data frame for survival analyses
df_base <- cluster_mapping %>%
  left_join(clinical, by = "bio_id_merge")
if ("FLT3_ITD_PCR.x" %in% names(df_base)) {
  df_base <- df_base %>%
    mutate(FLT3_ITD_PCR = coalesce(FLT3_ITD_PCR.x, FLT3_ITD_PCR.y)) %>%
    dplyr::select(-FLT3_ITD_PCR.x, -FLT3_ITD_PCR.y)
}

# DEGs (using adj.P <= 0.05 for robustness script, matching its original)
degs <- mitovsall %>%
  filter(logFC > 0, adj.P.Val <= 0.05, PG.Genes %in% rownames(vsn)) %>%
  pull(PG.Genes)

# HVP setup
protein_vars <- rowVars(vsn)
names(protein_vars) <- rownames(vsn)
vars_order <- order(protein_vars, decreasing = TRUE)
all_protein_names <- rownames(vsn)[vars_order]
n_q25 <- floor(nrow(vsn) * 0.25)
hvps_q25 <- all_protein_names[1:n_q25]

cat(sprintf("Gene set: %d DEGs ∩ %d HVPs = %d genes\n\n",
            length(degs), n_q25, length(mito_genes)))

# =============================================================================
# Helper functions
# =============================================================================

surv_p <- function(score_vec, df, cutoff_q = 0.75,
                   endpoint = "efs", subset = "all") {
  d <- df %>%
    mutate(score = score_vec[bio_id_merge]) %>%
    filter(!is.na(score))
  d$high <- d$score >= quantile(d$score, cutoff_q)

  if (subset == "flt3neg")      d <- d %>% filter(FLT3_ITD_PCR == 0)
  if (subset == "non_immature") d <- d %>% filter(cluster != "Immature_like")

  surv_col <- if (endpoint == "efs") "efs_days" else "os_days"
  stat_col <- if (endpoint == "efs") "efsstat"  else "stat"

  tryCatch({
    fit <- survdiff(as.formula(paste0("Surv(", surv_col, ",", stat_col, ") ~ high")),
                    data = d)
    1 - pchisq(fit$chisq, 1)
  }, error = function(e) NA)
}

cox_univar <- function(score_vec, df, cutoff_q = 0.75) {
  d <- df %>%
    mutate(score = score_vec[bio_id_merge],
           high  = score >= quantile(score, cutoff_q)) %>%
    filter(!is.na(score))
  tryCatch({
    s <- summary(coxph(Surv(efs_days, efsstat) ~ high, data = d))
    list(hr = s$conf.int["highTRUE", 1],
         ci_lo = s$conf.int["highTRUE", 3],
         ci_hi = s$conf.int["highTRUE", 4],
         p = s$coefficients["highTRUE", 5])
  }, error = function(e) list(hr = NA, ci_lo = NA, ci_hi = NA, p = NA))
}

cox_binary <- function(score_vec, df, cutoff_q = 0.75) {
  d <- df %>%
    mutate(score = score_vec[bio_id_merge],
           high  = score >= quantile(score, cutoff_q),
           cluster = factor(cluster)) %>%
    filter(!is.na(score))
  tryCatch({
    s <- summary(coxph(Surv(efs_days, efsstat) ~ high + FLT3_ITD_PCR + cluster, data = d))
    list(hr = s$conf.int["highTRUE", 1],
         ci_lo = s$conf.int["highTRUE", 3],
         ci_hi = s$conf.int["highTRUE", 4],
         p = s$coefficients["highTRUE", 5])
  }, error = function(e) list(hr = NA, ci_lo = NA, ci_hi = NA, p = NA))
}

cox_extended_flt3neg <- function(score_vec, df, cutoff_q = 0.75) {
  d <- df %>%
    mutate(score = score_vec[bio_id_merge],
           high  = score >= quantile(score, cutoff_q),
           cluster = factor(cluster),
           log_wbc = log10(pmax(wbc, 0.1)),
           log_ldh = log10(pmax(LDH_diag, 1)),
           gender = factor(gender),
           type_AML = factor(type_AML),
           treatment_ITT = factor(treatment_ITT)) %>%
    filter(FLT3_ITD_PCR == 0,
           !is.na(score), !is.na(log_wbc), !is.na(log_ldh),
           !is.na(age), !is.na(gender), !is.na(type_AML), !is.na(treatment_ITT))
  tryCatch({
    s <- summary(coxph(Surv(efs_days, efsstat) ~ high + cluster + log_wbc + log_ldh +
                         age + gender + type_AML + treatment_ITT, data = d))
    list(hr = s$conf.int["highTRUE", 1],
         ci_lo = s$conf.int["highTRUE", 3],
         ci_hi = s$conf.int["highTRUE", 4],
         p = s$coefficients["highTRUE", 5],
         n = nrow(d))
  }, error = function(e) list(hr = NA, ci_lo = NA, ci_hi = NA, p = NA, n = NA))
}

eval_full <- function(score_vec, df, cutoff_q = 0.75) {
  efs   <- surv_p(score_vec, df, cutoff_q, "efs", "all")
  flt3  <- surv_p(score_vec, df, cutoff_q, "efs", "flt3neg")
  os    <- surv_p(score_vec, df, cutoff_q, "os",  "all")
  ni    <- surv_p(score_vec, df, cutoff_q, "efs", "non_immature")
  f3d   <- cox_binary(score_vec, df, cutoff_q)
  fs3d  <- cox_extended_flt3neg(score_vec, df, cutoff_q)

  d <- df %>%
    mutate(score = score_vec[bio_id_merge]) %>%
    filter(!is.na(score))
  d$high <- d$score >= quantile(d$score, cutoff_q)
  ni_flt3_p <- tryCatch({
    fit <- survdiff(Surv(efs_days, efsstat) ~ high,
                    data = d %>% filter(cluster != "Immature_like" & FLT3_ITD_PCR == 0))
    1 - pchisq(fit$chisq, 1)
  }, error = function(e) NA)

  passes <- c(efs < 0.05, flt3 < 0.05, f3d$p < 0.05, fs3d$p < 0.05,
              ni < 0.05, ni_flt3_p < 0.05, os < 0.05)

  tibble(
    efs = efs, flt3neg = flt3, os = os,
    non_imm = ni, ni_flt3neg = ni_flt3_p,
    cox_hr = f3d$hr, cox_ci = sprintf("%.2f-%.2f", f3d$ci_lo, f3d$ci_hi),
    cox_p = f3d$p,
    ext_hr = fs3d$hr, ext_ci = sprintf("%.2f-%.2f", fs3d$ci_lo, fs3d$ci_hi),
    ext_p = fs3d$p,
    n_pass = sum(passes, na.rm = TRUE)
  )
}

# =============================================================================
# A. CUTOFF SENSITIVITY
# =============================================================================

cat("A. Cutoff sensitivity...\n")

cutoff_qs <- seq(0.50, 0.85, by = 0.05)
cutoff_res <- map_dfr(cutoff_qs, function(q) {
  d <- df_base %>% mutate(score = obs_score[bio_id_merge]) %>% filter(!is.na(score))
  n_high <- sum(d$score >= quantile(d$score, q))

  unadj <- cox_univar(obs_score, df_base, q)
  adj   <- cox_binary(obs_score, df_base, q)

  tibble(
    quantile = q, n_high = n_high,
    efs     = surv_p(obs_score, df_base, q, "efs", "all"),
    flt3neg = surv_p(obs_score, df_base, q, "efs", "flt3neg"),
    os      = surv_p(obs_score, df_base, q, "os",  "all"),
    non_imm = surv_p(obs_score, df_base, q, "efs", "non_immature"),
    unadj_hr = unadj$hr, unadj_ci_lo = unadj$ci_lo,
    unadj_ci_hi = unadj$ci_hi, unadj_p = unadj$p,
    cox_hr = adj$hr, cox_ci_lo = adj$ci_lo,
    cox_ci_hi = adj$ci_hi, cox_p = adj$p
  )
})

# =============================================================================
# B. GENE SET PERMUTATION TEST
# =============================================================================

cat("B. Permutation test (2,000 random gene sets)...\n")

set.seed(42)
n_perm <- 2000

obs_efs  <- surv_p(obs_score, df_base, 0.75, "efs", "all")
obs_flt3 <- surv_p(obs_score, df_base, 0.75, "efs", "flt3neg")
obs_os   <- surv_p(obs_score, df_base, 0.75, "os",  "all")

perm_pool <- setdiff(rownames(vsn), mito_genes)
perm_efs <- perm_flt3 <- perm_os <- numeric(n_perm)

for (i in seq_len(n_perm)) {
  random_genes <- sample(perm_pool, length(mito_genes))
  perm_ss <- gsva(ssgseaParam(vsn, list(random = random_genes)), verbose = FALSE)
  perm_score <- setNames(as.numeric(perm_ss[1, ]), colnames(perm_ss))

  perm_efs[i]  <- surv_p(perm_score, df_base, 0.75, "efs", "all")
  perm_flt3[i] <- surv_p(perm_score, df_base, 0.75, "efs", "flt3neg")
  perm_os[i]   <- surv_p(perm_score, df_base, 0.75, "os",  "all")

  if (i %% 100 == 0) cat(sprintf("  %d/%d\n", i, n_perm))
}

perm_p_efs  <- mean(perm_efs  <= obs_efs)
perm_p_flt3 <- mean(perm_flt3 <= obs_flt3)
perm_p_os   <- mean(perm_os   <= obs_os)

cat(sprintf("  Permutation p: EFS=%.3f, FLT3-=%.3f, OS=%.3f\n",
            perm_p_efs, perm_p_flt3, perm_p_os))

# =============================================================================
# C. BOOTSTRAP CONFIDENCE INTERVALS
# =============================================================================

cat("C. Bootstrap CIs (1,000 resamples)...\n")

set.seed(123)
n_boot <- 1000
boot_efs <- boot_flt3 <- boot_os <- numeric(n_boot)

for (i in seq_len(n_boot)) {
  boot_idx <- sample(nrow(df_base), replace = TRUE)
  boot_df  <- df_base[boot_idx, ]
  boot_df$score <- obs_score[boot_df$bio_id_merge]
  cutoff <- quantile(boot_df$score, 0.75, na.rm = TRUE)
  boot_df$high <- boot_df$score >= cutoff

  boot_efs[i] <- tryCatch({
    fit <- survdiff(Surv(efs_days, efsstat) ~ high, data = boot_df)
    1 - pchisq(fit$chisq, 1)
  }, error = function(e) NA)

  boot_flt3[i] <- tryCatch({
    fit <- survdiff(Surv(efs_days, efsstat) ~ high,
                    data = boot_df %>% filter(FLT3_ITD_PCR == 0))
    1 - pchisq(fit$chisq, 1)
  }, error = function(e) NA)

  boot_os[i] <- tryCatch({
    fit <- survdiff(Surv(os_days, stat) ~ high, data = boot_df)
    1 - pchisq(fit$chisq, 1)
  }, error = function(e) NA)
}

cat(sprintf("  Bootstrap p<0.05: EFS=%.1f%%, FLT3-=%.1f%%, OS=%.1f%%\n",
            mean(boot_efs < 0.05, na.rm = TRUE) * 100,
            mean(boot_flt3 < 0.05, na.rm = TRUE) * 100,
            mean(boot_os < 0.05, na.rm = TRUE) * 100))

# =============================================================================
# D. DEG FILTER SENSITIVITY
# =============================================================================

cat("D. DEG filter sensitivity...\n")

# D1: Varying adj.P (logFC > 0 fixed)
adjp_grid <- c(0.001, 0.01, 0.05, 0.1, 0.2)
deg_adjp_res <- map_dfr(adjp_grid, function(pval) {
  test_degs  <- mitovsall %>%
    filter(logFC > 0, adj.P.Val <= pval, PG.Genes %in% rownames(vsn)) %>%
    pull(PG.Genes)
  test_genes <- intersect(test_degs, hvps_q25)
  if (length(test_genes) < 10)
    return(tibble(adjp = pval, n_degs = length(test_degs),
                  n_genes = length(test_genes),
                  efs = NA, flt3neg = NA, os = NA, cox_p = NA))
  test_ss    <- gsva(ssgseaParam(vsn, list(mito = test_genes)), verbose = FALSE)
  test_score <- setNames(as.numeric(test_ss[1, ]), colnames(test_ss))
  tibble(
    adjp = pval, n_degs = length(test_degs), n_genes = length(test_genes),
    efs     = surv_p(test_score, df_base, 0.75, "efs", "all"),
    flt3neg = surv_p(test_score, df_base, 0.75, "efs", "flt3neg"),
    os      = surv_p(test_score, df_base, 0.75, "os",  "all"),
    cox_p   = cox_binary(test_score, df_base, 0.75)$p
  )
})

# D2: Varying logFC (adj.P <= 0.1 fixed)
lfc_grid <- c(0.0, 0.25, 0.5, 0.75, 1.0)
deg_lfc_res <- map_dfr(lfc_grid, function(lfc) {
  test_degs  <- mitovsall %>%
    filter(logFC > lfc, adj.P.Val <= 0.1, PG.Genes %in% rownames(vsn)) %>%
    pull(PG.Genes)
  test_genes <- intersect(test_degs, hvps_q25)
  if (length(test_genes) < 10)
    return(tibble(logFC = lfc, n_degs = length(test_degs),
                  n_genes = length(test_genes),
                  efs = NA, flt3neg = NA, os = NA, cox_p = NA))
  test_ss    <- gsva(ssgseaParam(vsn, list(mito = test_genes)), verbose = FALSE)
  test_score <- setNames(as.numeric(test_ss[1, ]), colnames(test_ss))
  tibble(
    logFC = lfc, n_degs = length(test_degs), n_genes = length(test_genes),
    efs     = surv_p(test_score, df_base, 0.75, "efs", "all"),
    flt3neg = surv_p(test_score, df_base, 0.75, "efs", "flt3neg"),
    os      = surv_p(test_score, df_base, 0.75, "os",  "all"),
    cox_p   = cox_binary(test_score, df_base, 0.75)$p
  )
})

deg_sensitivity <- bind_rows(
  deg_adjp_res %>% mutate(param = "adj.P", value = adjp) %>% dplyr::select(-adjp),
  deg_lfc_res  %>% mutate(param = "logFC", value = logFC) %>% dplyr::select(-logFC)
)

# =============================================================================
# E. HVP THRESHOLD / SCORING METHOD COMPARISON
# =============================================================================

cat("E. HVP threshold & scoring method comparison...\n")

target <- setNames(cluster_mapping$mito_score, cluster_mapping$bio_id_merge)

hvp_grid <- c(1068, 1100, 1200, 1300, 1500, n_q25, 2000)
hvp_grid <- sort(unique(hvp_grid))

hvp_res <- list()

for (hvp_n in hvp_grid) {
  hvps <- all_protein_names[1:hvp_n]
  genes <- intersect(degs, hvps)

  for (method in c("ssGSEA", "colMeans")) {
    label <- sprintf("%s_HVP%d", method, hvp_n)

    if (method == "ssGSEA") {
      ss <- gsva(ssgseaParam(vsn, list(mito = genes)), verbose = FALSE)
      score <- setNames(as.numeric(ss[1, ]), colnames(ss))
    } else {
      score <- setNames(colMeans2(vsn[genes, ]), colnames(vsn))
    }

    ev <- eval_full(score, df_base, 0.75)
    r  <- cor(score[names(target)], target[names(score)])

    hvp_res[[label]] <- ev %>%
      mutate(label = label, method = method, hvp_n = hvp_n,
             n_genes = length(genes), r_target = r)
  }
  cat(sprintf("  HVP=%d (n=%d genes): done\n", hvp_n, length(genes)))
}

# Special thresholds
special_thresholds <- list(
  "var_gt_mean"     = names(which(protein_vars > mean(protein_vars))),
  "var_gt_mean05sd" = names(which(protein_vars > mean(protein_vars) + 0.5 * sd(protein_vars))),
  "all_DEGs"        = rownames(vsn)
)

for (st_name in names(special_thresholds)) {
  genes <- intersect(degs, special_thresholds[[st_name]])
  label <- paste0("ssGSEA_", st_name)

  ss <- gsva(ssgseaParam(vsn, list(mito = genes)), verbose = FALSE)
  score <- setNames(as.numeric(ss[1, ]), colnames(ss))
  ev <- eval_full(score, df_base, 0.75)
  r  <- cor(score[names(target)], target[names(score)])

  hvp_res[[label]] <- ev %>%
    mutate(label = label, method = "ssGSEA",
           hvp_n = length(special_thresholds[[st_name]]),
           n_genes = length(genes), r_target = r)
  cat(sprintf("  %s (n=%d genes): done\n", st_name, length(genes)))
}

hvp_comparison <- bind_rows(hvp_res)

# =============================================================================
# F. CONTINUOUS COX REGRESSION
# =============================================================================

cat("F. Continuous Cox regression...\n")

df_cox <- df_base %>%
  mutate(score   = obs_score[bio_id_merge],
         score_z = as.numeric(scale(score)),
         cluster = factor(cluster),
         log_wbc = log10(pmax(wbc, 0.1)),
         log_ldh = log10(pmax(LDH_diag, 1)),
         gender  = factor(gender),
         type_AML = factor(type_AML),
         treatment_ITT = factor(treatment_ITT)) %>%
  filter(!is.na(score))

continuous_cox_res <- list()

for (ep in c("efs", "os")) {
  surv_col <- if (ep == "efs") "efs_days" else "os_days"
  stat_col <- if (ep == "efs") "efsstat"  else "stat"
  fit <- coxph(as.formula(paste0("Surv(", surv_col, ",", stat_col, ") ~ score_z")),
               data = df_cox)
  s <- summary(fit)
  cat(sprintf("  Univariate %s: HR=%.3f (%.3f-%.3f), p=%.4f\n", toupper(ep),
              s$conf.int["score_z", 1], s$conf.int["score_z", 3],
              s$conf.int["score_z", 4], s$coefficients["score_z", 5]))
  continuous_cox_res[[paste0("univar_", ep)]] <- list(
    hr = s$conf.int["score_z", 1], ci_lo = s$conf.int["score_z", 3],
    ci_hi = s$conf.int["score_z", 4], p = s$coefficients["score_z", 5])
}

# Multivariate
fit <- coxph(Surv(efs_days, efsstat) ~ score_z + FLT3_ITD_PCR + cluster, data = df_cox)
s <- summary(fit)
continuous_cox_res[["multivar_efs"]] <- list(
  hr = s$conf.int["score_z", 1], ci_lo = s$conf.int["score_z", 3],
  ci_hi = s$conf.int["score_z", 4], p = s$coefficients["score_z", 5])

# Extended FLT3-neg
df_f3 <- df_cox %>%
  filter(FLT3_ITD_PCR == 0, !is.na(log_wbc), !is.na(log_ldh),
         !is.na(age), !is.na(gender), !is.na(type_AML), !is.na(treatment_ITT))
fit <- coxph(Surv(efs_days, efsstat) ~ score_z + cluster + log_wbc + log_ldh +
               age + gender + type_AML + treatment_ITT, data = df_f3)
s <- summary(fit)
continuous_cox_res[["extended_flt3neg"]] <- list(
  hr = s$conf.int["score_z", 1], ci_lo = s$conf.int["score_z", 3],
  ci_hi = s$conf.int["score_z", 4], p = s$coefficients["score_z", 5])

# Spline nonlinearity test
tryCatch({
  library(rms)
  dd <- datadist(df_cox); options(datadist = "dd")
  fit_sp <- cph(Surv(efs_days, efsstat) ~ rcs(score_z, 3) + FLT3_ITD_PCR + cluster,
                data = df_cox, x = TRUE, y = TRUE)
  a <- anova(fit_sp)
  nl_row <- grep("[Nn]onlinear", rownames(a))
  if (length(nl_row) > 0) {
    cat(sprintf("  Spline nonlinearity p: %.4f\n", a[nl_row[1], "P"]))
    continuous_cox_res[["nonlinear_p"]] <- a[nl_row[1], "P"]
  }
}, error = function(e) cat(sprintf("  Spline test skipped: %s\n", e$message)))

# =============================================================================
# G. LEAVE-ONE-OUT PROTEIN STABILITY
# =============================================================================

cat("G. Leave-one-out stability...\n")

n_genes <- length(mito_genes)
loo_results <- tibble()

for (i in seq_along(mito_genes)) {
  remaining <- mito_genes[-i]
  loo_ss    <- gsva(ssgseaParam(vsn, list(mito = remaining)), verbose = FALSE)
  loo_score <- setNames(as.numeric(loo_ss[1, ]), colnames(loo_ss))

  loo_results <- bind_rows(loo_results, tibble(
    dropped_gene = mito_genes[i],
    r_with_full  = cor(loo_score, obs_score[names(loo_score)]),
    efs_p        = surv_p(loo_score, df_base, 0.75, "efs", "all"),
    flt3neg_p    = surv_p(loo_score, df_base, 0.75, "efs", "flt3neg"),
    os_p         = surv_p(loo_score, df_base, 0.75, "os",  "all"),
    cox_p        = cox_binary(loo_score, df_base, 0.75)$p
  ))

  if (i %% 50 == 0) cat(sprintf("  %d/%d\n", i, n_genes))
}

loo_results <- loo_results %>%
  mutate(all_sig = efs_p < 0.05 & flt3neg_p < 0.05 & os_p < 0.05 & cox_p < 0.05)

cat(sprintf("  All 4 criteria retained: %d/%d (%.1f%%)\n",
            sum(loo_results$all_sig), n_genes, 100 * mean(loo_results$all_sig)))

# =============================================================================
# H. SPLIT-HALF CROSS-VALIDATION
# =============================================================================

cat("H. Split-half cross-validation (500 splits)...\n")

set.seed(42)
n_splits  <- 500
half_size <- floor(n_genes / 2)
split_results <- tibble()

for (i in seq_len(n_splits)) {
  idx    <- sample(n_genes, half_size)
  half_a <- mito_genes[idx]
  half_b <- mito_genes[-idx]

  ss_a <- gsva(ssgseaParam(vsn, list(mito = half_a)), verbose = FALSE)
  ss_b <- gsva(ssgseaParam(vsn, list(mito = half_b)), verbose = FALSE)
  score_a <- setNames(as.numeric(ss_a[1, ]), colnames(ss_a))
  score_b <- setNames(as.numeric(ss_b[1, ]), colnames(ss_b))

  split_results <- bind_rows(split_results, tibble(
    split  = i,
    r_ab   = cor(score_a, score_b),
    efs_a  = surv_p(score_a, df_base, 0.75, "efs", "all"),
    efs_b  = surv_p(score_b, df_base, 0.75, "efs", "all"),
    flt3_a = surv_p(score_a, df_base, 0.75, "efs", "flt3neg"),
    flt3_b = surv_p(score_b, df_base, 0.75, "efs", "flt3neg"),
    os_a   = surv_p(score_a, df_base, 0.75, "os",  "all"),
    os_b   = surv_p(score_b, df_base, 0.75, "os",  "all"),
    cox_a  = cox_binary(score_a, df_base, 0.75)$p,
    cox_b  = cox_binary(score_b, df_base, 0.75)$p
  ))

  if (i %% 50 == 0) cat(sprintf("  %d/%d\n", i, n_splits))
}

split_results <- split_results %>%
  mutate(
    both_efs  = efs_a < 0.05 & efs_b < 0.05,
    both_flt3 = flt3_a < 0.05 & flt3_b < 0.05,
    both_os   = os_a < 0.05 & os_b < 0.05,
    both_cox  = cox_a < 0.05 & cox_b < 0.05,
    either_efs  = efs_a < 0.05 | efs_b < 0.05,
    either_flt3 = flt3_a < 0.05 | flt3_b < 0.05,
    either_os   = os_a < 0.05 | os_b < 0.05,
    either_cox  = cox_a < 0.05 | cox_b < 0.05
  )

cat(sprintf("  Inter-half r: median=%.3f\n", median(split_results$r_ab)))
cat(sprintf("  Both halves p<0.05: EFS=%.1f%%, FLT3-=%.1f%%\n",
            mean(split_results$both_efs) * 100,
            mean(split_results$both_flt3) * 100))

# =============================================================================
# SAVE RESULTS
# =============================================================================

cat("\nSaving results...\n")

write.csv(cutoff_res, file.path(table_dir, "robustness_cutoff_sensitivity.csv"),
          row.names = FALSE)
write.csv(deg_sensitivity, file.path(table_dir, "robustness_deg_sensitivity.csv"),
          row.names = FALSE)
write.csv(hvp_comparison, file.path(table_dir, "robustness_hvp_comparison.csv"),
          row.names = FALSE)
write.csv(loo_results, file.path(table_dir, "robustness_loo.csv"),
          row.names = FALSE)
write.csv(split_results, file.path(table_dir, "robustness_split_half.csv"),
          row.names = FALSE)

saveRDS(list(
  cutoff      = cutoff_res,
  permutation = list(obs = c(efs = obs_efs, flt3 = obs_flt3, os = obs_os),
                     perm_efs = perm_efs, perm_flt3 = perm_flt3, perm_os = perm_os,
                     perm_p = c(efs = perm_p_efs, flt3 = perm_p_flt3, os = perm_p_os)),
  bootstrap   = list(boot_efs = boot_efs, boot_flt3 = boot_flt3, boot_os = boot_os),
  deg_sensitivity = deg_sensitivity,
  hvp_comparison  = hvp_comparison,
  continuous_cox  = continuous_cox_res,
  loo         = loo_results,
  split_half  = split_results,
  gene_set    = mito_genes,
  score       = obs_score
), file.path(output_dir, "robustness_results.RDS"))

# Table S2 (Excel)
tryCatch({
  writexl::write_xlsx(list(
    cutoff_sensitivity = cutoff_res,
    deg_sensitivity = deg_sensitivity,
    hvp_comparison = hvp_comparison,
    loo = loo_results,
    split_half = split_results
  ), path = file.path(table_dir, "TableS2.xlsx"))
  cat("Saved: TableS2.xlsx\n")
}, error = function(e) cat("writexl not available; skipping Excel output\n"))

# =============================================================================
# FIGURES
# =============================================================================

# --- Figure S3B: Permutation test histograms ---
p_perm_efs <- tibble(value = perm_efs) %>%
  ggplot(aes(value)) +
  geom_histogram(col = "black", fill = "lightgrey", bins = 40) +
  cowplot::theme_cowplot() +
  geom_vline(xintercept = obs_efs, col = "red") +
  annotate("text", x = 1, y = Inf, vjust = 1.5,
           label = sprintf("EFS observed p = %.3f\nPermutation p = %.3f",
                           obs_efs, perm_p_efs), hjust = 1) +
  labs(x = "EFS log-rank p (random gene sets)", y = "Frequency")

p_perm_os <- tibble(value = perm_os) %>%
  ggplot(aes(value)) +
  geom_histogram(col = "black", fill = "lightgrey", bins = 40) +
  cowplot::theme_cowplot() +
  geom_vline(xintercept = obs_os, col = "red") +
  annotate("text", x = 1, y = Inf, vjust = 1.5,
           label = sprintf("OS observed p = %.3f\nPermutation p = %.3f",
                           obs_os, perm_p_os), hjust = 1) +
  labs(x = "OS log-rank p (random gene sets)", y = "Frequency")

pdf(file.path(fig_dir, "FigS3B_permutation.pdf"), width = 8, height = 4)
cowplot::plot_grid(p_perm_efs, p_perm_os)
dev.off()

# --- Figure S3C: Cutoff sensitivity forest plot ---
df_forest <- cutoff_res %>%
  dplyr::select(quantile, n_high,
                unadj_hr, unadj_ci_lo, unadj_ci_hi, unadj_p,
                cox_hr, cox_ci_lo, cox_ci_hi, cox_p) %>%
  pivot_longer(
    cols = c(starts_with("unadj_"), starts_with("cox_")),
    names_to = c("model", ".value"),
    names_pattern = "(unadj|cox)_(.*)"
  ) %>%
  mutate(
    model = ifelse(model == "unadj", "Unadjusted", "Adjusted"),
    model = factor(model, levels = c("Unadjusted", "Adjusted")),
    quantile_label = sprintf("%.0f%%", quantile * 100),
    quantile_label = fct_reorder(quantile_label, quantile),
    sig = p < 0.05,
    chosen = quantile == 0.75,
    annotation = sprintf("%.2f (%.2f\u2013%.2f) p=%.3f", hr, ci_lo, ci_hi, p)
  )

p_cutoff <- ggplot(df_forest, aes(x = hr, y = quantile_label, color = model, shape = sig)) +
  annotate("rect",
           ymin = which(levels(df_forest$quantile_label) == "75%") - 0.4,
           ymax = which(levels(df_forest$quantile_label) == "75%") + 0.4,
           xmin = -Inf, xmax = Inf, fill = "grey90", alpha = 0.5) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_linerange(aes(xmin = ci_lo, xmax = ci_hi),
                 linewidth = 0.7, position = position_dodge(width = 0.5)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
  geom_text(aes(x = ci_hi, label = annotation),
            position = position_dodge(width = 0.5),
            hjust = -0.05, size = 2.4, show.legend = FALSE) +
  scale_x_log10(limits = c(0.4, 12)) +
  scale_color_manual(values = c("Unadjusted" = "#4393C3", "Adjusted" = "#D6604D")) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "p < 0.05", "FALSE" = "n.s."),
                     name = NULL) +
  labs(x = "Hazard Ratio (95% CI)", y = "Quantile cutoff", color = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    axis.title.y = element_text(margin = margin(r = 10)),
    plot.margin = margin(10, 20, 10, 10)
  )

ggsave(file.path(fig_dir, "FigS3C_cutoff_sensitivity.pdf"), p_cutoff,
       width = 5, height = 4)

cat("Done. Robustness results and figures saved.\n")
