#!/usr/bin/env Rscript

# Build the four-panel main pilot figure from saved results and cohort metadata.
# Statistics are computed upstream; this script changes presentation only.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# Accept explicit paths so the figure can be reproduced outside an interactive session.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("usage: 09b_make_manuscript_figure.R <repo_root> <atlas_dir> <figure_dir>")
}
root <- normalizePath(args[1], mustWork = TRUE)
atlas_dir <- args[2]
figure_dir <- args[3]
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# Use a restrained, color-blind-compatible palette and one manuscript-wide theme.
palette <- c(navy = "#24476B", teal = "#16817A", gold = "#D39C3F", red = "#B64B4B",
             purple = "#6F5AA8", grey = "#7B8790")
theme_paper <- theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, color = "grey30"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.key.height = grid::unit(0.45, "cm"),
        strip.background = element_rect(fill = "grey94", color = NA),
        strip.text = element_text(face = "bold"))

# Read the inferential result table and derive approximate moderated 95% intervals.
effects <- fread(file.path(atlas_dir, "module_effects.csv"))
effects[, standard_error := fifelse(is.finite(t) & abs(t) > 0, abs(logFC / t), NA_real_)]
effects[, lower_95 := logFC - 1.96 * standard_error]
effects[, upper_95 := logFC + 1.96 * standard_error]
effects[, fdr_state := factor(ifelse(module_fdr < 0.05, "FDR < 0.05", "FDR >= 0.05"),
                              levels = c("FDR < 0.05", "FDR >= 0.05"))]

# Human-readable module labels are fixed here and shared by all panels.
module_labels <- c(
  nk_identity = "NK identity",
  fcgr3a_adcc_effector = "FCGR3A/Fc signaling",
  cytotoxicity = "Cytotoxicity",
  ifng_chemokine = "IFN-gamma/chemokines",
  endothelial_identity = "Endothelial identity",
  endothelial_activation_injury = "Endothelial activation",
  tissue_injury_repair = "Tissue injury/repair",
  negative_control_t_cell = "T-cell control",
  negative_control_monocyte = "Monocyte control"
)
effects[, module_label := module_labels[module]]

# Panel A: display actual biopsy/sample counts used by the study rather than nominal accession size.
g50084 <- fread(file.path(root, "data/processed/GSE50084/samples.csv"))[compartment == "biopsy"]
g36059 <- fread(file.path(root, "data/processed/GSE36059/samples.csv"))
g275824 <- fread(file.path(root, "data/processed/GSE275824/samples.csv"))
counts <- rbindlist(list(
  g50084[, .N, by = .(group = analysis_group)][, dataset := "GSE50084"],
  g36059[analysis_group %in% c("ABMR", "TCMR", "non-rejecting"),
         .N, by = .(group = analysis_group)][, dataset := "GSE36059"],
  g275824[, .N, by = .(group = paste(treatment, timepoint, sep = " / "))][, dataset := "GSE275824"]
))
p_a <- ggplot(counts, aes(dataset, N, fill = dataset)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = N), vjust = -0.35, size = 3) +
  scale_fill_manual(values = c(GSE50084 = palette["teal"], GSE36059 = palette["navy"],
                               GSE275824 = palette["gold"]), guide = "none") +
  facet_wrap(~group, scales = "free_x", nrow = 2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "A  Public-cohort architecture", subtitle = "Observed profiles entering each analysis",
       x = NULL, y = "Profiles") + theme_paper +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7))

# Panel B: compare discovery, independent validation, and ABMR-versus-TCMR specificity directly.
primary_modules <- c("nk_identity", "fcgr3a_adcc_effector", "cytotoxicity",
                     "ifng_chemokine", "endothelial_activation_injury")
replication_contrasts <- c(
  ABMR_DSApos_vs_DSAneg_nonrejection = "Discovery: DSA+ ABMR vs DSA- NR",
  ABMR_vs_nonrejection = "Validation: ABMR vs NR",
  ABMR_vs_TCMR = "Specificity: ABMR vs TCMR"
)
replication <- effects[module %in% primary_modules & contrast %in% names(replication_contrasts)]
replication[, comparison := factor(replication_contrasts[contrast], levels = unname(replication_contrasts))]
replication[, module_label := factor(module_label, levels = rev(module_labels[primary_modules]))]
p_b <- ggplot(replication, aes(logFC, module_label, color = comparison)) +
  geom_vline(xintercept = 0, color = "grey75", linewidth = 0.45) +
  geom_errorbar(aes(xmin = lower_95, xmax = upper_95), width = 0.13, orientation = "y",
                position = position_dodge(width = 0.55)) +
  geom_point(aes(shape = fdr_state), size = 2.7, position = position_dodge(width = 0.55)) +
  scale_color_manual(values = c(palette["teal"], palette["navy"], palette["red"])) +
  scale_shape_manual(values = c("FDR < 0.05" = 16, "FDR >= 0.05" = 1)) +
  labs(title = "B  Replication does not imply ABMR specificity",
       subtitle = "Positive effects favor ABMR; negative specificity effects favor TCMR",
       x = "Standardized module effect (95% CI)", y = NULL, color = "Comparison", shape = "Evidence") +
  theme_paper + theme(legend.position = "bottom")

# Panel C: expose the difficult DSA-positive/no-rejection comparison instead of hiding its null result.
dsa_contrasts <- c(
  ABMR_DSApos_vs_DSAneg_nonrejection = "DSA+ ABMR vs DSA- NR",
  ABMR_DSApos_vs_DSApos_nonrejection = "DSA+ ABMR vs DSA+ NR",
  DSApos_nonrejection_vs_DSAneg_nonrejection = "DSA+ NR vs DSA- NR"
)
dsa <- effects[module %in% primary_modules & contrast %in% names(dsa_contrasts)]
dsa[, comparison := factor(dsa_contrasts[contrast], levels = unname(dsa_contrasts))]
dsa[, module_label := factor(module_label, levels = rev(module_labels[primary_modules]))]
p_c <- ggplot(dsa, aes(logFC, module_label, color = comparison)) +
  geom_vline(xintercept = 0, color = "grey75", linewidth = 0.45) +
  geom_errorbar(aes(xmin = lower_95, xmax = upper_95), width = 0.13, orientation = "y",
                position = position_dodge(width = 0.55)) +
  geom_point(aes(shape = fdr_state), size = 2.7, position = position_dodge(width = 0.55)) +
  scale_color_manual(values = c(palette["red"], palette["purple"], palette["gold"])) +
  scale_shape_manual(values = c("FDR < 0.05" = 16, "FDR >= 0.05" = 1)) +
  labs(title = "C  DSA exposure is not equivalent to rejection",
       subtitle = "The DSA+ non-rejection contrast remains an explicit falsification test",
       x = "Standardized module effect (95% CI)", y = NULL, color = "Comparison", shape = "Evidence") +
  theme_paper + theme(legend.position = "bottom")

# Panel D: show patient-blocked placebo-adjusted treatment effects at both protocol time points.
treatment <- effects[module %in% primary_modules & grepl("treatment_difference", contrast)]
treatment[, week := factor(ifelse(grepl("W24$", contrast), "Week 24", "Week 52"),
                           levels = c("Week 24", "Week 52"))]
treatment[, module_label := factor(module_label, levels = rev(module_labels[primary_modules]))]
p_d <- ggplot(treatment, aes(logFC, module_label, color = week)) +
  geom_vline(xintercept = 0, color = "grey75", linewidth = 0.45) +
  geom_errorbar(aes(xmin = lower_95, xmax = upper_95), width = 0.13, orientation = "y",
                position = position_dodge(width = 0.5)) +
  geom_point(aes(shape = fdr_state), size = 2.8, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("Week 24" = palette["teal"], "Week 52" = palette["red"])) +
  scale_shape_manual(values = c("FDR < 0.05" = 16, "FDR >= 0.05" = 1)) +
  labs(title = "D  Randomized treatment-associated module changes",
       subtitle = "Patient-blocked placebo-adjusted difference-in-differences",
       x = "Treatment effect (95% CI)", y = NULL, color = "Time point", shape = "Evidence") +
  theme_paper + theme(legend.position = "bottom")

# Assemble a balanced 2-by-2 figure with a concise, claim-calibrated footer.
main_figure <- (p_a | p_b) / (p_c | p_d) +
  plot_annotation(
    title = "MOSAIC-Tx: reference-aware qualification of ABMR molecular signatures",
    subtitle = "Three public cohorts; within-cohort standardization; BH FDR within each declared contrast",
    caption = paste0("NR, non-rejection. Effects are associative. Bulk expression cannot assign cell lineage or prove ",
                     "antibody engagement, target causality, or drug efficacy."),
    theme = theme(plot.title = element_text(face = "bold", size = 15),
                  plot.subtitle = element_text(size = 10),
                  plot.caption = element_text(size = 8, color = "grey35"))
  )

# Export high-resolution raster and editable vector versions at journal-compatible proportions.
ggsave(file.path(figure_dir, "main_figure_multicohort.png"), main_figure,
       width = 15.5, height = 11.5, dpi = 300, bg = "white")
ggsave(file.path(figure_dir, "main_figure_multicohort.pdf"), main_figure,
       width = 15.5, height = 11.5, device = "pdf")

# Save every panel independently to support editorial rearrangement and supplementary figures.
ggsave(file.path(figure_dir, "main_panel_A.png"), plot = p_a, width = 7.5, height = 5.5, dpi = 300, bg = "white")
ggsave(file.path(figure_dir, "main_panel_B.png"), plot = p_b, width = 7.5, height = 5.5, dpi = 300, bg = "white")
ggsave(file.path(figure_dir, "main_panel_C.png"), plot = p_c, width = 7.5, height = 5.5, dpi = 300, bg = "white")
ggsave(file.path(figure_dir, "main_panel_D.png"), plot = p_d, width = 7.5, height = 5.5, dpi = 300, bg = "white")
