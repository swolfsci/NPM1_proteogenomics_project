# =============================================================================
# 05_go_relapse_by_cluster.R
#
# Effect of gemtuzumab ozogamicin (GO) on relapse by proteomic differentiation
# cluster (post-hoc, exploratory).
#
# In AMLSG 09-09 the only significant benefit of GO was a lower cumulative
# incidence of relapse (CIR; Döhner et al., Lancet Haematol 2023), so relapse is
# the endpoint here. Conventions follow the trial report: CIR in patients who
# reached CR/CRi, time from remission, death in CR as competing risk, and
# cause-specific Cox models stratified by the randomisation stratum (age 18-60
# vs >60).
#
# Generates: Figure S (relapse by cluster and arm; GMP-like prognosis by arm)
#
# Input:  output/preprocessed_data.RDS (from 00_preprocessing.R)
#         output/differentiation_data.RDS (from 01_differentiation.R)
# Output: output/figures/FigS_cluster_relapse_by_arm.pdf
#         output/tables/GO_relapse_by_cluster.csv
#         output/tables/CIR_2y_by_cluster_and_arm.csv
#         output/tables/cluster_prognosis_by_arm_GMP.csv
#
# Requires clinical_data.csv with treatment arm, relapse and survival endpoints;
# skipped with a message if those columns are absent.
# =============================================================================

library(tidyverse)
library(survival)
library(patchwork)

source("R/utils.R")

# =============================================================================
# 1. Load data
# =============================================================================

cat("Loading data...\n")
preprocessed <- readRDS(file.path(output_dir, "preprocessed_data.RDS"))
clinical <- preprocessed$clinical

diff_data <- readRDS(file.path(output_dir, "differentiation_data.RDS"))
cluster_mapping <- diff_data$cluster_mapping

required_cols <- c("treatment_ITT", "age", "efs_days", "efsstat", "os_days", "stat",
                   "cuminc", "rfs_days")
missing_cols <- setdiff(required_cols, colnames(clinical))

if (length(missing_cols) > 0) {
  message("05_go_relapse_by_cluster.R: clinical columns missing (",
          paste(missing_cols, collapse = ", "), "); skipping.")
} else {

  # ===========================================================================
  # 2. Analysis cohort
  # ===========================================================================

  cluster_levels <- c("Immature_like", "GMP_like", "Commited_like", "Intermediate")

  df <- cluster_mapping %>%
    dplyr::select(bio_id_merge, cluster) %>%
    inner_join(dplyr::select(clinical, bio_id_merge, all_of(required_cols)),
               by = "bio_id_merge") %>%
    filter(!is.na(treatment_ITT), !is.na(efsstat)) %>%
    mutate(
      cluster    = factor(recode(cluster, Committed_like = "Commited_like"),
                          levels = cluster_levels),
      arm        = factor(if_else(str_detect(treatment_ITT, "GO"), "GO + ATRA", "ATRA"),
                          levels = c("ATRA", "GO + ATRA")),
      age_group  = if_else(age > 60, ">60", "18-60"),
      gmp_like   = as.numeric(cluster == "GMP_like"),
      efs_months = efs_days / 30.44,
      os_months  = os_days / 30.44
    )

  # Relapse is defined only for patients in CR/CRi (cuminc is NA otherwise);
  # rfs_days counts from the date of remission.
  relapse_df <- df %>%
    filter(!is.na(cuminc)) %>%
    mutate(
      rfs_months = rfs_days / 30.44,
      event      = factor(cuminc, levels = 0:2,
                          labels = c("censored", "relapse", "death_in_CR")),
      relapse    = as.numeric(cuminc == 1)
    )

  cat(sprintf("Cohort: %d patients; %d in CR/CRi with %d relapses and %d deaths in CR\n",
              nrow(df), nrow(relapse_df), sum(relapse_df$cuminc == 1),
              sum(relapse_df$cuminc == 2)))

  # --- Helpers ---

  hr_row <- function(fit, term) {
    s <- summary(fit)
    i <- which(rownames(s$coefficients) == term)
    tibble(HR = s$conf.int[i, 1], lower = s$conf.int[i, 3], upper = s$conf.int[i, 4],
           p = s$coefficients[i, "Pr(>|z|)"])
  }

  # Aalen-Johansen cumulative incidence of relapse (%) at 24 months
  cir_2y <- function(d) {
    fit <- survfit(Surv(rfs_months, event) ~ 1, data = d)
    s <- summary(fit, times = 24, extend = TRUE)
    100 * s$pstate[, which(fit$states == "relapse")]
  }

  lrt_p <- function(formula_reduced, formula_full, data) {
    anova(coxph(formula_reduced, data), coxph(formula_full, data))$`Pr(>|Chi|)`[2]
  }

  # ===========================================================================
  # 3. GO vs ATRA on relapse, overall and within each cluster
  # ===========================================================================

  cat("GO effect on relapse (cause-specific Cox, age-group stratified)...\n")

  overall <- hr_row(coxph(Surv(rfs_months, relapse) ~ arm + strata(age_group), relapse_df),
                    "armGO + ATRA")
  cat(sprintf("  all CR/CRi: HR %.2f [%.2f-%.2f], p = %.3f\n",
              overall$HR, overall$lower, overall$upper, overall$p))

  go_by_cluster <- relapse_df %>%
    group_by(cluster) %>%
    group_modify(~ bind_cols(
      tibble(n = nrow(.x),
             relapses_ATRA = sum(.x$relapse[.x$arm == "ATRA"]),
             relapses_GO   = sum(.x$relapse[.x$arm == "GO + ATRA"])),
      hr_row(coxph(Surv(rfs_months, relapse) ~ arm + strata(age_group), .x), "armGO + ATRA")
    )) %>%
    ungroup() %>%
    mutate(p_interaction_cluster_x_arm = lrt_p(
      Surv(rfs_months, relapse) ~ cluster + arm + strata(age_group),
      Surv(rfs_months, relapse) ~ cluster * arm + strata(age_group),
      relapse_df))

  print(go_by_cluster)
  write_csv(go_by_cluster, file.path(table_dir, "GO_relapse_by_cluster.csv"))

  cir_by_cluster_arm <- relapse_df %>%
    group_by(cluster, arm) %>%
    group_modify(~ tibble(n = nrow(.x), relapses = sum(.x$relapse), CIR_2y = cir_2y(.x))) %>%
    ungroup()

  write_csv(cir_by_cluster_arm, file.path(table_dir, "CIR_2y_by_cluster_and_arm.csv"))

  # ===========================================================================
  # 4. Prognostic effect of the GMP-like phenotype by treatment arm
  # ===========================================================================
  # The same treatment-by-cluster interaction, read from the prognostic side:
  # is the GMP-like advantage present in both arms?

  cat("GMP-like prognostic effect by arm...\n")

  endpoints <- list(
    EFS     = list(data = df,         surv = "Surv(efs_months, efsstat)"),
    OS      = list(data = df,         surv = "Surv(os_months, stat)"),
    Relapse = list(data = relapse_df, surv = "Surv(rfs_months, relapse)")
  )

  gmp_by_arm <- imap_dfr(endpoints, function(ep, name) {
    f <- function(rhs) as.formula(paste(ep$surv, "~", rhs))
    within <- map_dfr(c("Pooled", levels(df$arm)), function(a) {
      d <- if (a == "Pooled") ep$data else filter(ep$data, arm == a)
      hr_row(coxph(f("gmp_like + strata(age_group)"), d), "gmp_like") %>%
        mutate(arm = a, .before = 1)
    })
    int <- summary(coxph(f("gmp_like * arm + strata(age_group)"), ep$data))$coefficients
    within %>%
      mutate(
        endpoint = name, .before = 1,
        p_interaction_gmp_x_arm = int[grep("gmp_like:arm", rownames(int)), "Pr(>|z|)"],
        p_interaction_cluster_x_arm = lrt_p(f("cluster + arm + strata(age_group)"),
                                            f("cluster * arm + strata(age_group)"), ep$data)
      )
  })

  print(gmp_by_arm)
  write_csv(gmp_by_arm, file.path(table_dir, "cluster_prognosis_by_arm_GMP.csv"))

  # ===========================================================================
  # 5. Figure: relapse by cluster in each arm; GMP-like prognosis by arm
  # ===========================================================================

  cif_df <- map_dfr(levels(df$arm), function(a) {
    fit <- survfit(Surv(rfs_months, event) ~ cluster, data = filter(relapse_df, arm == a))
    k <- which(fit$states == "relapse")
    bind_rows(
      tibble(arm = a, cluster = sub("cluster=", "", names(fit$strata)), time = 0, cif = 0),
      tibble(arm = a, cluster = sub("cluster=", "", rep(names(fit$strata), fit$strata)),
             time = fit$time, cif = fit$pstate[, k])
    )
  }) %>%
    mutate(cluster = factor(cluster, levels = cluster_levels),
           arm = factor(arm, levels = levels(df$arm)))

  cif_labels <- cir_by_cluster_arm %>%
    arrange(arm, cluster) %>%
    group_by(arm) %>%
    mutate(y = 0.66 - 0.045 * (row_number() - 1),
           label = sprintf("%s  n=%d  2-y %.0f%%", cluster_labels[as.character(cluster)],
                           n, CIR_2y)) %>%
    ungroup()

  p_cif <- ggplot(cif_df, aes(time, cif, colour = cluster)) +
    geom_step(linewidth = 1.05) +
    geom_text(data = cif_labels, aes(x = 1, y = y, label = label, colour = cluster),
              hjust = 0, size = 3.1, fontface = "bold", show.legend = FALSE) +
    facet_wrap(~ arm) +
    scale_colour_manual(values = cluster_colors, labels = cluster_labels, name = NULL) +
    scale_y_continuous(labels = scales::percent) +
    coord_cartesian(xlim = c(0, 60), ylim = c(0, 0.68)) +
    labs(x = "Months from remission", y = "Cumulative incidence of relapse",
         title = "Relapse by differentiation cluster under standard therapy and with gemtuzumab ozogamicin",
         subtitle = sprintf("Patients in CR/CRi (n = %d); death in CR as competing risk. Curves shown to 60 months.",
                            nrow(relapse_df))) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "top",
          strip.background = element_rect(fill = "grey93", colour = NA),
          strip.text = element_text(face = "bold", size = 11),
          plot.title = element_text(face = "bold", size = 11.5),
          plot.subtitle = element_text(size = 8.8, colour = "grey35"))

  row_order <- paste(rep(c("EFS", "OS", "Relapse"), each = 3),
                     rep(c("Pooled", levels(df$arm)), 3), sep = "  |  ")
  interaction_text <- gmp_by_arm %>%
    distinct(endpoint, p_interaction_gmp_x_arm) %>%
    mutate(s = sprintf("%s %.3f", endpoint, p_interaction_gmp_x_arm)) %>%
    pull(s) %>%
    paste(collapse = " | ")

  p_forest <- gmp_by_arm %>%
    mutate(row = factor(paste(endpoint, arm, sep = "  |  "), levels = rev(row_order))) %>%
    ggplot(aes(HR, row, colour = arm)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
    geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.2, linewidth = 0.7,
                  orientation = "y") +
    geom_point(size = 2.8) +
    geom_text(aes(x = 4.2, label = sprintf("%.2f [%.2f-%.2f]  p=%.3f", HR, lower, upper, p)),
              hjust = 0, size = 3, colour = "grey20") +
    scale_colour_manual(values = c("Pooled" = "grey30", "ATRA" = "#2166ac",
                                   "GO + ATRA" = "#d6604d"), name = NULL) +
    scale_x_log10(limits = c(0.1, 30), breaks = c(0.25, 0.5, 1, 2, 4)) +
    labs(x = "Hazard ratio, GMP-like vs other clusters (log scale; < 1 = GMP-like better)",
         y = NULL,
         title = "Prognostic effect of the GMP-like phenotype by treatment arm",
         subtitle = paste0("Age-group stratified Cox models. GMP-like x arm interaction p: ",
                           interaction_text, ".")) +
    theme_bw(base_size = 10.5) +
    theme(panel.grid.minor = element_blank(), legend.position = "top",
          plot.title = element_text(face = "bold", size = 11.5),
          plot.subtitle = element_text(size = 8.8, colour = "grey35"))

  p_combined <- (p_cif / p_forest) +
    plot_layout(heights = c(1.25, 1)) +
    plot_annotation(
      tag_levels = "A",
      caption = paste("Post-hoc exploratory analysis. Differences in significance between arms",
                      "are not tests of interaction; see interaction p-values."),
      theme = theme(plot.caption = element_text(size = 8.5, colour = "grey35", hjust = 0))) &
    theme(plot.tag = element_text(face = "bold", size = 13))

  # Multi-panel PDFs written straight into a cloud-synced folder can end up as
  # empty files, so render to tempdir() and copy into place.
  tmp_pdf <- file.path(tempdir(), "FigS_cluster_relapse_by_arm.pdf")
  ggsave(tmp_pdf, p_combined, width = 11, height = 11, bg = "white")
  file.copy(tmp_pdf, file.path(fig_dir, "FigS_cluster_relapse_by_arm.pdf"), overwrite = TRUE)

  cat("Done. Saved: output/figures/FigS_cluster_relapse_by_arm.pdf\n")
}
