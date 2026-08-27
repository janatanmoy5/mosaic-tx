#!/usr/bin/env Rscript

# Render biopsy-level and patient-level results from the three downloaded cohorts.
# This script never simulates values: every point is linked to a GEO sample accession.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(jsonlite)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: 15_make_dataset_result_figures.R <repo_root> <output_dir>")
root <- normalizePath(args[1], mustWork = TRUE)
out_dir <- args[2]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cols <- c(control = "#707B83", difficult = "#D39C3F", abmr = "#B64B4B",
          tcmr = "#6F5AA8", placebo = "#24476B", treated = "#16817A")
theme_data <- theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9.5, color = "grey30"),
        strip.background = element_rect(fill = "grey94", color = NA),
        strip.text = element_text(face = "bold", size = 9),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom")

module_labels <- c(
  nk_identity = "NK identity", fcgr3a_adcc_effector = "FCGR3A/Fc signaling",
  cytotoxicity = "Cytotoxicity", ifng_chemokine = "IFN-gamma/chemokines",
  endothelial_identity = "Endothelial identity",
  endothelial_activation_injury = "Endothelial activation",
  tissue_injury_repair = "Tissue injury/repair",
  negative_control_t_cell = "T-cell control",
  negative_control_monocyte = "Monocyte control"
)
selected_modules <- c("fcgr3a_adcc_effector", "cytotoxicity", "ifng_chemokine",
                      "endothelial_activation_injury")

# Read one gene-by-sample matrix and enforce unique, uppercase gene symbols.
read_expression <- function(dataset) {
  x <- fread(file.path(root, "data/processed", dataset, "gene_expression.csv.gz"))
  setnames(x, 1L, "gene_symbol")
  x[, gene_symbol := toupper(gene_symbol)]
  x <- x[nzchar(gene_symbol)]
  if (anyDuplicated(x$gene_symbol)) {
    x <- x[, lapply(.SD, mean, na.rm = TRUE), by = gene_symbol]
  }
  x
}

# Calculate the prespecified mean gene-z module score used by the inference pipeline.
score_modules <- function(expr, sample_ids, dataset, metadata) {
  config <- fromJSON(file.path(root, "configs/mechanism_modules.json"), simplifyVector = FALSE)
  genes <- expr$gene_symbol
  mat <- as.matrix(expr[, ..sample_ids])
  storage.mode(mat) <- "double"
  rows <- lapply(names(config$modules), function(module) {
    declared <- toupper(unlist(config$modules[[module]]$genes))
    idx <- which(genes %chin% declared)
    if (length(idx) / length(declared) < config$rules$minimum_gene_coverage_fraction) return(NULL)
    block <- mat[idx, , drop = FALSE]
    z <- t(scale(t(block)))
    score <- colMeans(z, na.rm = TRUE)
    data.table(dataset = dataset, sample_accession = sample_ids,
               module = module, score = as.numeric(score))
  })
  ans <- rbindlist(rows)
  ans[metadata, on = "sample_accession"]
}

# PCA is computed only from observed expression, using the 1,000 most variable genes.
make_pca <- function(expr, sample_ids, metadata, dataset, group_col, group_levels,
                     group_labels, group_colors, title, subtitle) {
  mat <- as.matrix(expr[, ..sample_ids])
  storage.mode(mat) <- "double"
  variance <- apply(mat, 1L, var, na.rm = TRUE)
  keep <- head(order(variance, decreasing = TRUE, na.last = NA), min(1000L, sum(is.finite(variance))))
  fit <- prcomp(t(mat[keep, , drop = FALSE]), center = TRUE, scale. = TRUE)
  pc <- data.table(sample_accession = rownames(fit$x), PC1 = fit$x[, 1], PC2 = fit$x[, 2])
  pc <- metadata[pc, on = "sample_accession"]
  pc[, plot_group := factor(get(group_col), levels = group_levels, labels = group_labels)]
  fwrite(pc, file.path(out_dir, paste0(dataset, "_pca_coordinates.csv")))
  explained <- 100 * fit$sdev^2 / sum(fit$sdev^2)
  ggplot(pc, aes(PC1, PC2, color = plot_group)) +
    geom_hline(yintercept = 0, color = "grey90", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey90", linewidth = 0.3) +
    geom_point(size = 2.2, alpha = 0.82) +
    stat_ellipse(aes(group = plot_group), linewidth = 0.65, linetype = 2,
                 level = 0.80, show.legend = FALSE) +
    scale_color_manual(values = setNames(group_colors, group_labels), drop = FALSE) +
    labs(title = title, subtitle = subtitle,
         x = sprintf("PC1 (%.1f%%)", explained[1]), y = sprintf("PC2 (%.1f%%)", explained[2]),
         color = NULL) + theme_data
}

# Plot all observed sample scores with box summaries; no inferential claim is attached here.
make_distributions <- function(scores, group_col, group_levels, group_labels, group_colors,
                               title, subtitle) {
  x <- scores[module %in% selected_modules]
  x[, module_label := factor(module_labels[module], levels = module_labels[selected_modules])]
  x[, plot_group := factor(get(group_col), levels = group_levels, labels = group_labels)]
  ggplot(x, aes(plot_group, score, color = plot_group)) +
    geom_hline(yintercept = 0, color = "grey82", linewidth = 0.35) +
    geom_boxplot(width = 0.60, outlier.shape = NA, linewidth = 0.55) +
    geom_point(position = position_jitter(width = 0.13, height = 0, seed = 50084),
               alpha = 0.58, size = 1.35) +
    facet_wrap(~module_label, ncol = 2, scales = "free_y") +
    scale_color_manual(values = setNames(group_colors, group_labels), drop = FALSE) +
    labs(title = title, subtitle = subtitle, x = NULL,
         y = "Prespecified module score\n(mean within-cohort gene z-score)", color = NULL) +
    theme_data + theme(axis.text.x = element_text(angle = 18, hjust = 1), legend.position = "none")
}

# GSE50084: restrict to the 61 primary allograft biopsies; blood samples are excluded.
m50084 <- fread(file.path(root, "data/processed/GSE50084/samples.csv"))[compartment == "biopsy"]
e50084 <- read_expression("GSE50084")
s50084 <- score_modules(e50084, m50084$sample_accession, "GSE50084", m50084)
g50084_levels <- c("DSA_negative_control", "DSA_positive_no_rejection", "DSA_positive_rejection")
g50084_labels <- c("DSA- / no rejection\ncontrol (n=20)", "DSA+ / no rejection\ncontrol (n=13)",
                   "DSA+ / rejection\nABMR-enriched (n=28)")
p_a <- make_pca(e50084, m50084$sample_accession, m50084, "GSE50084", "analysis_group",
                g50084_levels, g50084_labels, unname(cols[c("control", "difficult", "abmr")]),
                "A  GSE50084: 61 primary biopsies",
                "Unsupervised PCA of 1,000 most-variable genes; color labels are metadata-derived")
p_b <- make_distributions(s50084, "analysis_group", g50084_levels, g50084_labels,
                          unname(cols[c("control", "difficult", "abmr")]),
                          "B  GSE50084: DSA state is separated from rejection",
                          "Each point is one primary biopsy; DSA+ non-rejection is the difficult control")

# GSE36059: retain the predeclared non-rejecting, TCMR, and ABMR specificity groups.
m36059 <- fread(file.path(root, "data/processed/GSE36059/samples.csv"))[
  analysis_group %in% c("non-rejecting", "TCMR", "ABMR")]
e36059 <- read_expression("GSE36059")
s36059 <- score_modules(e36059, m36059$sample_accession, "GSE36059", m36059)
g36059_levels <- c("non-rejecting", "TCMR", "ABMR")
g36059_labels <- c("Non-rejecting\ncontrol (n=281)", "TCMR\ndisease control (n=35)", "ABMR\n(n=65)")
p_c <- make_pca(e36059, m36059$sample_accession, m36059, "GSE36059", "analysis_group",
                g36059_levels, g36059_labels, unname(cols[c("control", "tcmr", "abmr")]),
                "C  GSE36059: independent replication and specificity",
                "Unsupervised PCA; TCMR is retained as an explicit disease-control class")
p_d <- make_distributions(s36059, "analysis_group", g36059_levels, g36059_labels,
                          unname(cols[c("control", "tcmr", "abmr")]),
                          "D  GSE36059: observed biopsy-level module scores",
                          "All 381 prespecified validation biopsies are displayed")

# GSE275824: visualize actual longitudinal patient measurements and randomized assignment.
m275824 <- fread(file.path(root, "data/processed/GSE275824/samples.csv"))
score_labels <- c(abmract = "ABMR activity", abmrifng = "ABMR IFN-gamma",
                  abmrnk = "ABMR NK", abmrendo = "ABMR endothelium")
trajectory <- melt(m275824,
                   id.vars = c("sample_accession", "patient_id", "treatment", "timepoint"),
                   measure.vars = names(score_labels), variable.name = "published_score",
                   value.name = "score")
trajectory[, published_score := factor(score_labels[as.character(published_score)],
                                       levels = unname(score_labels))]
trajectory[, timepoint := factor(timepoint, levels = c("Baseline", "Week24", "Week52"))]
trajectory[, treatment_label := factor(treatment, levels = c("Placebo", "Felzartamab"),
                                       labels = c("Placebo (control)", "Felzartamab (treated)"))]
fwrite(trajectory, file.path(out_dir, "GSE275824_patient_trajectories.csv"))
mean_trajectory <- trajectory[, .(mean = mean(score, na.rm = TRUE),
                                  se = sd(score, na.rm = TRUE) / sqrt(sum(is.finite(score)))),
                              by = .(published_score, treatment_label, timepoint)]
p_e <- ggplot(trajectory, aes(timepoint, score, group = interaction(patient_id, treatment_label),
                              color = treatment_label)) +
  geom_line(alpha = 0.22, linewidth = 0.48) +
  geom_point(alpha = 0.38, size = 1.25) +
  geom_errorbar(data = mean_trajectory,
                aes(x = timepoint, y = mean, ymin = mean - 1.96 * se,
                    ymax = mean + 1.96 * se, color = treatment_label,
                    group = treatment_label), inherit.aes = FALSE,
                width = 0.07, linewidth = 0.72, alpha = 0.95) +
  geom_line(data = mean_trajectory,
            aes(x = timepoint, y = mean, color = treatment_label, group = treatment_label),
            inherit.aes = FALSE, linewidth = 1.15) +
  geom_point(data = mean_trajectory,
             aes(x = timepoint, y = mean, color = treatment_label),
             inherit.aes = FALSE, size = 2.6) +
  facet_wrap(~published_score, nrow = 1, scales = "free_y") +
  scale_color_manual(values = c("Placebo (control)" = unname(cols["placebo"]),
                                "Felzartamab (treated)" = unname(cols["treated"]))) +
  labs(title = "E  GSE275824: randomized longitudinal response",
       subtitle = "63 biopsies from 22 patients; thin lines are patients, thick lines are group means ± approximate 95% CI",
       x = NULL, y = "Published molecular score", color = "Randomized arm") + theme_data

# Export audit tables containing every expression-derived point used in panels B and D.
audit_scores <- rbindlist(list(s50084, s36059), fill = TRUE)
fwrite(audit_scores, file.path(out_dir, "individual_sample_module_scores.csv"))

# Assemble one sample-level verification figure; do not export redundant standalone panels.
main <- (p_a | p_b) / (p_c | p_d) / p_e +
  plot_layout(heights = c(1, 1, 0.88)) +
  plot_annotation(
    title = "MOSAIC-Tx: observed multicohort evidence for ABMR-associated molecular programs",
    subtitle = "No simulated data: 61 GSE50084 biopsies, 381 prespecified GSE36059 biopsies, and 63 GSE275824 longitudinal biopsies",
    caption = paste0("PCA and distributions are descriptive. Group labels come from public metadata. ",
                     "Bulk RNA cannot assign a program to a specific cell lineage or establish causality."),
    theme = theme(plot.title = element_text(face = "bold", size = 16),
                  plot.subtitle = element_text(size = 10),
                  plot.caption = element_text(size = 8.5, color = "grey30")))
ggsave(file.path(out_dir, "dataset_results_main.png"), main, width = 16, height = 19,
       dpi = 300, bg = "white")
ggsave(file.path(out_dir, "dataset_results_main.pdf"), main, width = 16, height = 19,
       device = "pdf")

cat("Dataset-level figures generated from", nrow(m50084), "GSE50084 biopsies,",
    nrow(m36059), "GSE36059 biopsies, and", nrow(m275824), "GSE275824 biopsies.\n")
