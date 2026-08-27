#!/usr/bin/env Rscript

# Purpose: render publication-style pilot figures only from saved, auditable result tables.
# The script does not recompute statistics; this separates inference from presentation.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# Require the atlas directory and figure directory as explicit arguments.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: 09_make_pilot_figures.R <atlas_dir> <figure_dir>")
atlas_dir <- args[1]
figure_dir <- args[2]
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# Load the complete module effect table and define readable plotting labels.
effects <- fread(file.path(atlas_dir, "module_effects.csv"))
contrast_labels <- c(
  ABMR_DSApos_vs_DSAneg_nonrejection = "DSA+ ABMR vs DSA- non-rejection",
  ABMR_DSApos_vs_DSApos_nonrejection = "DSA+ ABMR vs DSA+ non-rejection",
  DSApos_nonrejection_vs_DSAneg_nonrejection = "DSA+ non-rejection vs DSA- non-rejection",
  ABMR_vs_nonrejection = "ABMR vs non-rejection (validation)",
  ABMR_vs_TCMR = "ABMR vs TCMR (specificity)",
  treatment_difference_in_differences_W24 = "Felzartamab DeltaDelta at week 24",
  treatment_difference_in_differences_W52 = "Felzartamab DeltaDelta at week 52"
)
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
effects[, contrast_label := factor(contrast_labels[contrast], levels = rev(unname(contrast_labels)))]
effects[, module_label := factor(module_labels[module], levels = unname(module_labels))]
effects[, significance := fifelse(module_fdr < 0.05, "FDR < 0.05", "FDR >= 0.05")]

# Panel A: the complete cross-cohort effect-size atlas; dot size communicates adjusted evidence.
p_atlas <- ggplot(effects, aes(contrast_label, module_label, color = logFC,
                               size = pmin(-log10(pmax(module_fdr, 1e-12)), 6))) +
  geom_point() +
  scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                        name = "Module effect") +
  scale_size_continuous(range = c(2, 8), name = "-log10 FDR\n(capped at 6)") +
  coord_flip() +
  labs(title = "A  Cross-cohort ABMR mechanism-module atlas", x = NULL, y = NULL,
       subtitle = "Scores are standardized within cohort; positive values favor the first group") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1),
        plot.title = element_text(face = "bold"))

# Panel B: focus on the randomized longitudinal immune-versus-tissue module trajectories.
treatment <- effects[grepl("treatment_difference", contrast)]
treatment[, week := factor(ifelse(grepl("W24$", contrast), "Week 24", "Week 52"),
                           levels = c("Week 24", "Week 52"))]
p_treatment <- ggplot(treatment, aes(logFC, module_label, color = week)) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.5) +
  geom_errorbar(aes(xmin = logFC - 1.96 * abs(logFC / t), xmax = logFC + 1.96 * abs(logFC / t)),
                width = 0.15, orientation = "y", position = position_dodge(width = 0.5)) +
  geom_point(aes(shape = significance), size = 3, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("Week 24" = "#007C83", "Week 52" = "#A44A3F")) +
  labs(title = "B  Placebo-adjusted longitudinal module effects", x = "Treatment difference-in-differences",
       y = NULL, color = NULL, shape = NULL,
       caption = "Bars are approximate 95% confidence intervals from limma moderated standard errors.") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

# Combine panels into a single journal-width figure and save both vector and raster formats.
combined <- p_atlas / p_treatment + plot_layout(heights = c(1.35, 1))
ggsave(file.path(figure_dir, "pilot_signature_robustness_atlas.pdf"), combined,
       width = 14, height = 12, device = "pdf")
ggsave(file.path(figure_dir, "pilot_signature_robustness_atlas.png"), combined,
       width = 14, height = 12, dpi = 300, bg = "white")

# Export individual panels for manuscript-layout flexibility.
ggsave(file.path(figure_dir, "pilot_panel_a_cross_cohort_atlas.png"), p_atlas,
       width = 13, height = 7, dpi = 300, bg = "white")
ggsave(file.path(figure_dir, "pilot_panel_b_treatment_effects.png"), p_treatment,
       width = 10, height = 6, dpi = 300, bg = "white")
