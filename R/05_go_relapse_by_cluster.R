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
# Beyond the discrete clusters, the GO effect is tested along the continuous
# differentiation scores and diffusion components, and protein by protein
# (protein x arm interaction screen with pathway-level GSEA).
#
# Generates: Figure S (relapse by cluster and arm; GMP-like prognosis by arm)
#            Figure S (GO effect on relapse along DC1/DC2; per-decile estimates)
#
# Input:  output/preprocessed_data.RDS (from 00_preprocessing.R)
#         output/differentiation_data.RDS (from 01_differentiation.R)
# Output: output/figures/FigS_cluster_relapse_by_arm.pdf
#         output/figures/FigS_GO_relapse_HR_curve_DC.pdf
#         output/figures/FigS_GO_relapse_HR_by_DC_deciles.pdf
#         output/tables/GO_relapse_by_cluster.csv
#         output/tables/CIR_2y_by_cluster_and_arm.csv
#         output/tables/cluster_prognosis_by_arm_GMP.csv
#         output/tables/GO_relapse_interaction_differentiation.csv
#         output/tables/GO_relapse_HR_curve_DC_tests.csv
#         output/tables/GO_relapse_HR_by_DC_deciles.csv
#         output/tables/GO_CIR_protein_interaction_screen.csv
#         output/tables/GO_CIR_protein_interaction_GSEA.csv
#
# Requires clinical_data.csv with treatment arm, relapse and survival endpoints;
# skipped with a message if those columns are absent. The decile estimates need
# coxphf and the pathway analysis needs fgsea and msigdbr; each is skipped with
# a message if its package is missing.
# =============================================================================

library(tidyverse)
library(survival)
library(splines)
library(patchwork)

source("R/utils.R")

# =============================================================================
# 1. Load data
# =============================================================================

cat("Loading data...\n")
preprocessed <- readRDS(file.path(output_dir, "preprocessed_data.RDS"))
clinical <- preprocessed$clinical

vsn <- preprocessed$vsn_matrix

diff_data <- readRDS(file.path(output_dir, "differentiation_data.RDS"))
cluster_mapping <- diff_data$cluster_mapping
diff_scores_wide <- diff_data$diff_scores_wide

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
    dplyr::select(bio_id_merge, cluster, DC1, DC2) %>%
    left_join(diff_scores_wide, by = "bio_id_merge") %>%
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

  save_pdf <- function(p, name, width, height) {
    tmp <- file.path(tempdir(), name)
    ggsave(tmp, p, width = width, height = height, bg = "white")
    file.copy(tmp, file.path(fig_dir, name), overwrite = TRUE)
  }

  # ===========================================================================
  # 6. GO effect along the continuous differentiation axes
  # ===========================================================================
  # Each continuous score and diffusion component is tested for an interaction
  # with treatment arm: a straight-line interaction (HR per SD) and a natural
  # spline (3 df) interaction that allows any shape along the axis (LRT).

  cat("GO x continuous differentiation interaction...\n")

  surv_rel <- "Surv(rfs_months, relapse)"

  add_spline <- function(d, x) {
    basis <- ns(x, df = 3)
    bind_cols(d, as_tibble(unclass(basis)[, 1:3], .name_repair = ~ paste0("b", 1:3)))
  }

  test_axis <- function(v) {
    d <- relapse_df %>%
      filter(!is.na(.data[[v]])) %>%
      mutate(z = scale(.data[[v]])[, 1])
    d <- add_spline(d, d$z)
    sl <- summary(coxph(as.formula(paste(surv_rel, "~ arm * z + strata(age_group)")), d))$coef
    tibble(
      variable             = v,
      HR_interaction_perSD = exp(sl[grep(":z$", rownames(sl)), "coef"]),
      p_linear             = sl[grep(":z$", rownames(sl)), "Pr(>|z|)"],
      p_spline_any         = lrt_p(
        as.formula(paste(surv_rel, "~ arm + b1 + b2 + b3 + strata(age_group)")),
        as.formula(paste(surv_rel, "~ arm * (b1 + b2 + b3) + strata(age_group)")), d)
    )
  }

  go_interaction <- map_dfr(c("Immature_like", "Committed_like", "GMP_like", "DC1", "DC2"),
                            test_axis) %>%
    bind_rows(tibble(variable = "Cluster (4 levels)",
                     p_linear = unique(go_by_cluster$p_interaction_cluster_x_arm)))

  print(go_interaction)
  write_csv(go_interaction, file.path(table_dir, "GO_relapse_interaction_differentiation.csv"))

  # --- HR curves along DC1 and DC2 (spline interaction model) ---

  curves <- list()
  curve_tests <- list()
  for (ax in c("DC1", "DC2")) {
    d <- relapse_df %>% mutate(z = .data[[ax]])
    basis <- ns(d$z, df = 3)
    d <- add_spline(d, d$z)
    f_add <- coxph(Surv(rfs_months, relapse) ~ arm + b1 + b2 + b3 + strata(age_group), d)
    f_lin <- coxph(Surv(rfs_months, relapse) ~ arm * z + b1 + b2 + b3 + strata(age_group), d)
    f_spl <- coxph(Surv(rfs_months, relapse) ~ arm * (b1 + b2 + b3) + strata(age_group), d)
    curve_tests[[ax]] <- tibble(
      axis                    = ax,
      p_any_interaction       = anova(f_add, f_spl)$`Pr(>|Chi|)`[2],
      p_linear_interaction    = anova(f_add, f_lin)$`Pr(>|Chi|)`[2],
      p_departure_from_linear = anova(f_lin, f_spl)$`Pr(>|Chi|)`[2]
    )

    # log-HR of GO vs ATRA at each grid point = arm coefficient + arm:spline terms
    grid <- seq(quantile(d$z, 0.03), quantile(d$z, 0.97), length.out = 150)
    b_grid <- predict(basis, grid)
    nm <- names(coef(f_spl))
    C <- matrix(0, length(grid), length(nm), dimnames = list(NULL, nm))
    C[, grep("^arm", nm)[1]] <- 1
    for (j in 1:3) C[, grep(paste0(":b", j, "$"), nm)] <- b_grid[, j]
    lhr <- drop(C %*% coef(f_spl))
    se <- sqrt(rowSums((C %*% vcov(f_spl)) * C))
    curves[[ax]] <- tibble(axis = ax, x = grid, HR = exp(lhr),
                           lo = exp(lhr - 1.96 * se), hi = exp(lhr + 1.96 * se))
  }
  curve_tests <- bind_rows(curve_tests)
  curves <- bind_rows(curves)

  print(curve_tests)
  write_csv(curve_tests, file.path(table_dir, "GO_relapse_HR_curve_DC_tests.csv"))

  # --- Per-decile estimates (Firth-penalised; few events per decile) ---

  if (requireNamespace("coxphf", quietly = TRUE)) {
    deciles <- map_dfr(c("DC1", "DC2"), function(ax) {
      relapse_df %>%
        mutate(decile = ntile(.data[[ax]], 10)) %>%
        group_by(decile) %>%
        group_modify(function(d, k) {
          f <- coxphf::coxphf(Surv(rfs_months, relapse) ~ arm, data = as.data.frame(d))
          tibble(n = nrow(d),
                 relapses_ATRA = sum(d$relapse[d$arm == "ATRA"]),
                 relapses_GO   = sum(d$relapse[d$arm == "GO + ATRA"]),
                 HR = exp(unname(f$coefficients[1])),
                 lower = unname(f$ci.lower[1]), upper = unname(f$ci.upper[1]),
                 p = unname(f$prob[1]),
                 cluster = names(which.max(table(d$cluster))),
                 x = median(d[[ax]]))
        }) %>%
        ungroup() %>%
        mutate(axis = ax, .before = 1)
    })
    write_csv(deciles, file.path(table_dir, "GO_relapse_HR_by_DC_deciles.csv"))
  } else {
    message("coxphf not installed; skipping per-decile GO estimates.")
    deciles <- NULL
  }

  # --- Figures ---

  y_lim <- c(0.03, 30)
  hr_breaks <- c(0.05, 0.1, 0.25, 0.5, 1, 2, 4, 10, 25)
  strip <- setNames(
    sprintf("%s: effect varies along axis p = %.3f; straight-line trend p = %.3f; curvature p = %.3f",
            curve_tests$axis, curve_tests$p_any_interaction,
            curve_tests$p_linear_interaction, curve_tests$p_departure_from_linear),
    curve_tests$axis)
  facet_lab <- function(a) factor(strip[a], levels = strip)
  hr_theme <- theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "top",
          strip.background = element_rect(fill = "grey93", colour = NA),
          strip.text = element_text(face = "bold", hjust = 0, size = 9.5),
          plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 8.8, colour = "grey35"),
          plot.caption = element_text(size = 8.3, colour = "grey35", hjust = 0))
  y_lab <- "Cause-specific HR for relapse, GO + ATRA vs ATRA (log scale)"

  rug <- bind_rows(transmute(relapse_df, axis = "DC1", x = DC1),
                   transmute(relapse_df, axis = "DC2", x = DC2))

  p_curve <- ggplot() +
    geom_hline(yintercept = 1, colour = "grey45") +
    geom_ribbon(data = curves %>% mutate(f = facet_lab(axis), lo = pmax(lo, y_lim[1]),
                                         hi = pmin(hi, y_lim[2])),
                aes(x, ymin = lo, ymax = hi), fill = "#d6604d", alpha = 0.15) +
    geom_line(data = curves %>% mutate(f = facet_lab(axis)), aes(x, HR),
              colour = "#d6604d", linewidth = 1.1) +
    geom_rug(data = rug %>% mutate(f = facet_lab(axis)), aes(x = x),
             alpha = 0.35, length = unit(0.02, "npc")) +
    facet_wrap(~ f, ncol = 1, scales = "free_x") +
    scale_y_log10(breaks = hr_breaks) +
    coord_cartesian(ylim = y_lim) +
    labs(x = "Diffusion component value (left = low, right = high)", y = y_lab,
         title = "Effect of gemtuzumab ozogamicin on relapse along DC1 and DC2",
         subtitle = sprintf(paste0("Patients in CR/CRi (n = %d, %d relapses). Line and band: HR and 95%% CI from an age-stratified\n",
                                   "cause-specific Cox model with a natural spline (3 df) for each axis and its interaction with arm."),
                            nrow(relapse_df), sum(relapse_df$relapse)),
         caption = "Post-hoc exploratory analysis. Tests are likelihood-ratio tests of the arm interaction terms.") +
    hr_theme

  if (!is.null(deciles)) {
    p_curve <- p_curve +
      geom_point(data = deciles %>%
                   mutate(f = facet_lab(axis), HR = pmin(pmax(HR, y_lim[1]), y_lim[2]),
                          cluster = factor(cluster, levels = names(cluster_colors))),
                 aes(x, HR, fill = cluster), shape = 21, size = 3, colour = "grey15") +
      scale_fill_manual(values = cluster_colors, labels = cluster_labels, drop = FALSE,
                        name = "Decile estimate (Firth), dominant cluster")

    p_deciles <- deciles %>%
      mutate(f = facet_lab(axis), lower = pmax(lower, y_lim[1]), upper = pmin(upper, y_lim[2]),
             cluster = factor(cluster, levels = names(cluster_colors)),
             txt = sprintf("%d/%d", relapses_ATRA, relapses_GO)) %>%
      ggplot(aes(decile, HR)) +
      geom_hline(yintercept = 1, colour = "grey45") +
      geom_hline(yintercept = overall$HR, colour = "#d6604d", linetype = 2) +
      geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25, colour = "grey40") +
      geom_point(aes(fill = cluster), shape = 21, size = 3.6, colour = "grey15") +
      geom_text(aes(y = y_lim[1] * 1.25, label = txt), size = 2.8, colour = "grey30") +
      facet_wrap(~ f, ncol = 1) +
      scale_fill_manual(values = cluster_colors, labels = cluster_labels, drop = FALSE,
                        name = "Dominant cluster in decile") +
      scale_y_log10(breaks = hr_breaks) +
      scale_x_continuous(breaks = 1:10, labels = c("1\nlowest", 2:9, "10\nhighest")) +
      coord_cartesian(ylim = y_lim) +
      labs(x = "Decile", y = y_lab,
           title = "Effect of gemtuzumab ozogamicin on relapse across deciles of DC1 and DC2",
           subtitle = sprintf(paste0("Patients in CR/CRi (n = %d, %d relapses). Dashed line: overall age-stratified HR %.2f.\n",
                                     "Numbers above the axis are relapses in the ATRA / GO arm. Error bars are 95%% CIs, clipped at the panel edge."),
                              nrow(relapse_df), sum(relapse_df$relapse), overall$HR),
           caption = paste0("Per-decile estimates from Firth-penalised cause-specific Cox models with arm as the only covariate.\n",
                            "Single-decile estimates rest on few events and are unreliable on their own.")) +
      hr_theme

    save_pdf(p_deciles, "FigS_GO_relapse_HR_by_DC_deciles.pdf", 9.5, 9)
  }

  save_pdf(p_curve, "FigS_GO_relapse_HR_curve_DC.pdf", 9.5, 9)

  # ===========================================================================
  # 7. Protein x GO interaction screen on relapse
  # ===========================================================================
  # One age-stratified cause-specific Cox model per protein (z-scored across
  # the CR/CRi cohort) with a protein x arm interaction. HR_int < 1 means that
  # patients with higher protein levels derive a larger relapse benefit from GO.

  cat("Protein x GO interaction screen...\n")

  cr_ids <- intersect(relapse_df$bio_id_merge, colnames(vsn))
  cr <- relapse_df[match(cr_ids, relapse_df$bio_id_merge), ]
  X <- t(scale(t(vsn[, cr_ids])))

  screen <- map_dfr(seq_len(nrow(X)), function(i) {
    x <- X[i, ]
    f <- tryCatch(coxph(Surv(cr$rfs_months, cr$relapse) ~ x * cr$arm + strata(cr$age_group)),
                  error = function(e) NULL)
    if (is.null(f)) return(tibble(protein = rownames(X)[i], z = NA_real_, p = NA_real_,
                                  HR_int = NA_real_))
    sc <- summary(f)$coef
    r <- grep("^x:", rownames(sc))
    tibble(protein = rownames(X)[i], z = sc[r, "z"], p = sc[r, "Pr(>|z|)"],
           HR_int = exp(sc[r, "coef"]))
  }) %>%
    mutate(q = p.adjust(p, "BH")) %>%
    arrange(p)

  screen_ok <- filter(screen, !is.na(p))
  cat(sprintf("  %d proteins; FDR < 0.05: %d; FDR < 0.20: %d; raw p < 0.05: %d (expected %.0f)\n",
              nrow(screen_ok), sum(screen_ok$q < 0.05), sum(screen_ok$q < 0.20),
              sum(screen_ok$p < 0.05), 0.05 * nrow(screen_ok)))
  cat(sprintf("  genomic-control lambda = %.2f\n",
              median(screen_ok$z^2) / qchisq(0.5, 1)))
  cat("  A priori candidates (CD33 = GO target, BCL2, MCL1):\n")
  print(filter(screen_ok, protein %in% c("CD33", "BCL2", "MCL1")))

  write_csv(screen, file.path(table_dir, "GO_CIR_protein_interaction_screen.csv"))

  # --- Pathway level: GSEA on the interaction z statistics ---

  if (requireNamespace("fgsea", quietly = TRUE) && requireNamespace("msigdbr", quietly = TRUE)) {
    msig <- bind_rows(
      msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
      msigdbr::msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME"))
    set.seed(1)
    go_gsea <- fgsea::fgsea(split(msig$gene_symbol, msig$gs_name),
                            sort(setNames(screen_ok$z, screen_ok$protein)),
                            minSize = 10, maxSize = 500) %>%
      as_tibble() %>%
      arrange(padj) %>%
      mutate(leadingEdge = map_chr(leadingEdge, ~ paste(head(.x, 20), collapse = "/")))
    cat(sprintf("  GSEA on interaction z: %d pathways with padj < 0.05\n",
                sum(go_gsea$padj < 0.05, na.rm = TRUE)))
    write_csv(go_gsea, file.path(table_dir, "GO_CIR_protein_interaction_GSEA.csv"))
  } else {
    message("fgsea/msigdbr not installed; skipping GSEA of the interaction screen.")
  }

  cat("Done. Saved GO relapse figures to output/figures/ and tables to output/tables/\n")
}
