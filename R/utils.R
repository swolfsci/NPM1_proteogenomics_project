# =============================================================================
# Shared helper functions
# =============================================================================

library(tidyverse)
library(survival)
library(survminer)

# --- Paths -------------------------------------------------------------------

data_dir   <- "data"
output_dir <- "output"
fig_dir    <- file.path(output_dir, "figures")
table_dir  <- file.path(output_dir, "tables")

dir.create(fig_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# --- Cluster colors and labels -----------------------------------------------

cluster_colors <- c(
  "Commited_like"  = "#219ebc",
  "GMP_like"       = "#ffb703",
  "Immature_like"  = "#023047",
  "Intermediate"   = "#5bb14fff"
)

cluster_labels <- c(
  "Commited_like"  = "Committed-like",
  "GMP_like"       = "GMP-like",
  "Immature_like"  = "Immature-like",
  "Intermediate"   = "Intermediate"
)

# --- Differentiation score aggregation ---------------------------------------

aggregate_differentiation_scores <- function(gsva_long) {
  gsva_long %>%
    mutate(
      signature_group = case_when(
        name %in% c("LSPC-Primed-Top100", "MLL_LSC_Somervaille2009_Up") ~ "Immature_like",
        name %in% c("furtwaengler_gmdp") ~ "GMP_like",
        TRUE ~ "Committed_like"
      )
    ) %>%
    summarize(value = mean(value), .by = c(bio_id_merge, signature_group)) %>%
    pivot_wider(id_cols = bio_id_merge, names_from = signature_group, values_from = value)
}

# --- Genotype / phenotype association ----------------------------------------

# Wilcoxon rank-sum test of each binary mutation against the continuous
# differentiation scores, with BH correction across the gene x score family.
# Used for the final Fig 2D panel (02_mutations.R) and the genetic-vs-non-genetic
# decomposition (05_reviewer_analyses.R).
test_genotype_differentiation <- function(diff_scores_wide, clinical, mutation_cols,
                                          score_cols = c("Committed_like",
                                                         "Immature_like", "GMP_like")) {
  diff_scores_wide %>%
    left_join(clinical, by = "bio_id_merge") %>%
    dplyr::select(all_of(score_cols), any_of(mutation_cols)) %>%
    pivot_longer(cols = any_of(mutation_cols),
                 names_to = "gene", values_to = "mutation_status") %>%
    pivot_longer(cols = all_of(score_cols),
                 names_to = "name", values_to = "value") %>%
    filter(!is.na(mutation_status)) %>%
    group_by(gene, name) %>%
    rstatix::wilcox_test(value ~ mutation_status) %>%
    rstatix::adjust_pvalue(method = "BH") %>%
    rstatix::add_significance()
}

# --- Bootstrap R^2 for variance decomposition --------------------------------

# Row-resampling bootstrap of the R^2 of lm(outcome ~ .) on a data frame that
# contains `outcome` plus a set of predictor columns. Returns the boot object.
bootstrap_r2 <- function(data, outcome, R = 2000, seed = 42) {
  set.seed(seed)
  stat <- function(d, idx) {
    f <- as.formula(paste(outcome, "~ ."))
    summary(lm(f, data = d[idx, ]))$r.squared
  }
  boot::boot(data, statistic = stat, R = R)
}

# --- GSEA helpers (used in 03_mito_score.R) ----------------------------------

run_limma_contrast <- function(expr, annot, group_col, contrast_str) {
  design <- model.matrix(~ 0 + factor(annot[[group_col]]))
  colnames(design) <- levels(factor(annot[[group_col]]))
  fit <- limma::lmFit(expr, design) %>%
    limma::contrasts.fit(limma::makeContrasts(contrasts = contrast_str, levels = design)) %>%
    limma::eBayes()
  limma::topTable(fit, number = Inf, sort.by = "none") %>%
    rownames_to_column("protein") %>%
    as_tibble()
}

run_gsea_all <- function(dep_df, hallmark_df, mitocarta_df, zeng_df) {
  ranked <- dep_df %>%
    dplyr::select(protein, t) %>%
    deframe() %>%
    sort(decreasing = TRUE)

  gsea_h <- clusterProfiler::GSEA(ranked, TERM2GENE = hallmark_df,
                                   pvalueCutoff = 1, minGSSize = 10, maxGSSize = 500,
                                   eps = 0, verbose = FALSE)
  gsea_go <- clusterProfiler::gseGO(ranked, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                     keyType = "SYMBOL", ont = "BP",
                                     pvalueCutoff = 1, minGSSize = 15,
                                     maxGSSize = 500, eps = 0, verbose = FALSE)
  gsea_mc <- clusterProfiler::GSEA(ranked, TERM2GENE = mitocarta_df,
                                    pvalueCutoff = 1, minGSSize = 5, maxGSSize = 500,
                                    eps = 0, verbose = FALSE)
  gsea_z <- clusterProfiler::GSEA(ranked, TERM2GENE = zeng_df,
                                   pvalueCutoff = 1, minGSSize = 5, maxGSSize = 500,
                                   eps = 0, verbose = FALSE)

  bind_rows(
    gsea_h@result  %>% as_tibble() %>% mutate(source = "Hallmark"),
    gsea_go@result %>% as_tibble() %>% mutate(source = "GO:BP"),
    gsea_mc@result %>% as_tibble() %>% mutate(source = "MitoCarta"),
    gsea_z@result  %>% as_tibble() %>% mutate(source = "Zeng")
  )
}

dedup_by_jaccard <- function(df, cutoff = 0.5) {
  if (nrow(df) <= 1) return(df)
  genes <- str_split(df$core_enrichment, "/")
  keep <- rep(TRUE, nrow(df))
  for (i in 2:nrow(df)) {
    for (j in 1:(i - 1)) {
      if (!keep[j]) next
      jac <- length(intersect(genes[[i]], genes[[j]])) /
        length(union(genes[[i]], genes[[j]]))
      if (jac > cutoff) { keep[i] <- FALSE; break }
    }
  }
  df[keep, ]
}

select_paired_pathways <- function(gsea_c1, gsea_c2, n_per_source = 5) {
  wide <- gsea_c1 %>%
    select(ID, Description, source, NES_c1 = NES, padj_c1 = p.adjust,
           core_enrichment_c1 = core_enrichment) %>%
    full_join(
      gsea_c2 %>% select(ID, Description, source, NES_c2 = NES, padj_c2 = p.adjust,
                         core_enrichment_c2 = core_enrichment),
      by = c("ID", "source")
    ) %>%
    mutate(Description = coalesce(Description.x, Description.y)) %>%
    select(-Description.x, -Description.y) %>%
    filter(padj_c1 < 0.05 | padj_c2 < 0.05) %>%
    mutate(max_abs_NES = pmax(abs(NES_c1), abs(NES_c2), na.rm = TRUE))

  go_candidates <- wide %>%
    filter(source == "GO:BP") %>%
    arrange(-max_abs_NES) %>%
    head(n_per_source * 3)

  if (nrow(go_candidates) > 1) {
    genes <- map2(
      str_split(replace_na(go_candidates$core_enrichment_c1, ""), "/"),
      str_split(replace_na(go_candidates$core_enrichment_c2, ""), "/"),
      ~ union(.x[.x != ""], .y[.y != ""])
    )
    keep <- rep(TRUE, nrow(go_candidates))
    for (i in 2:nrow(go_candidates)) {
      for (j in 1:(i - 1)) {
        if (!keep[j]) next
        jac <- length(intersect(genes[[i]], genes[[j]])) /
          length(union(genes[[i]], genes[[j]]))
        if (jac > 0.5) { keep[i] <- FALSE; break }
      }
    }
    go_candidates <- go_candidates[keep, ]
  }
  go_selected <- go_candidates %>% head(n_per_source)

  other_selected <- wide %>%
    filter(source != "GO:BP") %>%
    group_by(source) %>%
    slice_max(max_abs_NES, n = n_per_source, with_ties = FALSE) %>%
    ungroup()

  bind_rows(other_selected, go_selected)
}

select_top_pathways <- function(gsea_df, n_per_source = 5) {
  sig <- gsea_df %>% filter(p.adjust < 0.05)

  go_sig <- sig %>%
    filter(source == "GO:BP") %>%
    slice_max(abs(NES), n = n_per_source * 3, with_ties = FALSE) %>%
    dedup_by_jaccard(cutoff = 0.5) %>%
    head(n_per_source)

  other_sig <- sig %>%
    filter(source != "GO:BP") %>%
    group_by(source) %>%
    slice_max(abs(NES), n = n_per_source, with_ties = FALSE) %>%
    ungroup()

  bind_rows(other_sig, go_sig)
}

clean_labels <- function(df) {
  df %>%
    mutate(
      label = case_when(
        source == "Hallmark" ~ str_remove(ID, "^HALLMARK_") %>%
          str_replace_all("_", " ") %>% str_to_title() %>%
          str_replace("\\bE2f\\b", "E2F") %>%
          str_replace("\\bMyc\\b", "MYC") %>%
          str_replace("\\bDna\\b", "DNA") %>%
          str_replace("\\bUv\\b", "UV") %>%
          str_replace("\\bG2m\\b", "G2M") %>%
          str_replace("\\bTnfa\\b", "TNFa") %>%
          str_replace("\\bIfn\\b", "IFN") %>%
          str_replace("\\bIl(\\d)", "IL\\1") %>%
          str_replace("\\bP53\\b", "p53") %>%
          str_replace("\\bKras\\b", "KRAS") %>%
          str_replace("\\bMtorc1\\b", "mTORC1") %>%
          str_replace("\\bTgf\\b", "TGF") %>%
          str_replace("\\bNfkb\\b", "NF-kB") %>%
          str_replace("\\bWnt\\b", "WNT") %>%
          str_replace("\\bRos\\b", "ROS") %>%
          str_replace("\\bPi3k\\b", "PI3K") %>%
          str_replace("\\bAkt\\b", "AKT"),
        source == "GO:BP" ~ str_to_sentence(Description),
        source == "MitoCarta" ~ ID,
        source == "Zeng" ~ str_remove(ID, "-Top100$"),
        TRUE ~ ID
      ),
      label = str_trunc(label, 45)
    )
}

# --- Plot functions ----------------------------------------------------------

make_paired_dotplot <- function(selected_wide, gsea_c1, gsea_c2,
                                 title, subtitle,
                                 c1_label = "vs Immature high-mito",
                                 c2_label = "vs Non-immature reference") {

  contrast_colors <- c("#E64B35", "#4DBBD5")
  names(contrast_colors) <- c(c1_label, c2_label)

  source_labels <- c(
    "Hallmark"  = "MSigDB Hallmark",
    "GO:BP"     = "GO Biological Process",
    "MitoCarta" = "MitoCarta 3.0",
    "Zeng"      = "Cell type signatures (Zeng)"
  )

  plot_data <- selected_wide %>%
    clean_labels() %>%
    select(ID, source, label, NES_c1, padj_c1, NES_c2, padj_c2) %>%
    pivot_longer(
      cols = c(NES_c1, NES_c2),
      names_to = "contrast_var",
      values_to = "NES"
    ) %>%
    mutate(
      padj = ifelse(contrast_var == "NES_c1", padj_c1, padj_c2),
      contrast = ifelse(contrast_var == "NES_c1", c1_label, c2_label),
      NES = replace_na(NES, 0),
      padj = replace_na(padj, 1),
      significant = padj < 0.05
    ) %>%
    select(ID, source, label, contrast, NES, padj, significant)

  pathway_summary <- selected_wide %>%
    clean_labels() %>%
    mutate(max_nes = pmax(abs(NES_c1), abs(NES_c2), na.rm = TRUE))

  source_order <- c("Hallmark", "MitoCarta", "GO:BP", "Zeng")
  pathway_summary <- pathway_summary %>%
    mutate(source = factor(source, levels = source_order)) %>%
    arrange(source, -max_nes)

  label_order <- pathway_summary$label
  plot_data <- plot_data %>%
    mutate(
      label = factor(label, levels = rev(label_order)),
      contrast = factor(contrast, levels = c(c1_label, c2_label)),
      source = factor(source, levels = source_order)
    )

  ggplot(plot_data, aes(x = NES, y = label)) +
    geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4) +
    geom_line(aes(group = label), color = "grey70", linewidth = 0.5) +
    geom_point(aes(color = contrast, shape = significant, size = significant)) +
    facet_grid(source ~ ., scales = "free_y", space = "free_y", switch = "y",
               labeller = labeller(source = source_labels)) +
    scale_color_manual(values = contrast_colors, name = "Contrast") +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                       labels = c("TRUE" = "adj. p < 0.05", "FALSE" = "n.s."),
                       name = "Significance") +
    scale_size_manual(values = c("TRUE" = 3, "FALSE" = 2), guide = "none") +
    labs(x = "Normalized Enrichment Score", y = NULL,
         title = title, subtitle = subtitle) +
    theme_bw(base_size = 10) +
    theme(
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", size = 9),
      strip.background = element_rect(fill = "grey95", color = NA),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal",
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 8, color = "grey40")
    )
}

make_lollipop_plot <- function(top_df, title, subtitle,
                                up_label = "Enriched in high-mito",
                                down_label = "Depleted in high-mito") {
  source_colors <- c(
    "Hallmark"  = "#D6604D",
    "GO:BP"     = "#4393C3",
    "MitoCarta" = "#74ADD1",
    "Zeng"      = "#8DA0CB"
  )

  plot_df <- top_df %>%
    mutate(
      direction = ifelse(NES > 0, up_label, down_label),
      direction = factor(direction, levels = c(up_label, down_label))
    ) %>%
    arrange(direction, NES) %>%
    mutate(label = fct_inorder(label))

  ggplot(plot_df, aes(x = NES, y = label)) +
    geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4) +
    geom_segment(aes(x = 0, xend = NES, yend = label, color = source),
                 linewidth = 0.7) +
    geom_point(aes(color = source, size = -log10(p.adjust))) +
    facet_grid(direction ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_color_manual(values = source_colors, name = "Collection") +
    scale_size_continuous(range = c(2, 5), name = expression(-log[10](adj.~italic(p)))) +
    labs(x = "Normalized Enrichment Score", y = NULL,
         title = title, subtitle = subtitle) +
    theme_bw(base_size = 10) +
    theme(
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold", size = 9),
      strip.background = element_rect(fill = "grey95", color = NA),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 8, color = "grey40")
    )
}
