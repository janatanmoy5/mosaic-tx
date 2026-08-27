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
theme_data <- theme_classic(base_size=9.5) +
  theme(plot.title=element_text(face="bold",size=10.5,margin=margin(b=2)),
        plot.subtitle=element_text(size=8,color="grey30",margin=margin(b=2)),
        strip.background=element_rect(fill="grey94",color=NA),
        strip.text=element_text(face="bold",size=8.5),
        axis.title=element_text(face="bold",size=9.5),
        axis.title.x=element_text(margin=margin(t=2)),
        axis.title.y=element_text(margin=margin(r=2)),
        axis.text=element_text(size=8.5,color="black"),
        axis.text.y=element_text(margin=margin(r=1)),
        legend.position="bottom",
        legend.title=element_text(face="bold",size=7.5),
        legend.text=element_text(size=7),
        legend.key.height=grid::unit(8,"pt"),
        legend.key.width=grid::unit(10,"pt"),
        legend.spacing.x=grid::unit(2,"pt"),
        legend.background=element_rect(fill="white",color="grey45",linewidth=.3),
        legend.box.background=element_rect(fill="white",color="grey45",linewidth=.3),
        legend.margin=margin(2,3,2,3),
        plot.margin=margin(3,3,3,3))

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

# GSE50084 includes biopsy and whole blood. Primary inference is biopsy-only.
m50084_all <- fread(file.path(root, "data/processed/GSE50084/samples.csv"))
m50084 <- m50084_all[compartment == "biopsy"]
e50084 <- read_expression("GSE50084")
s50084 <- score_modules(e50084, m50084$sample_accession, "GSE50084", m50084)
g50084_levels <- c("DSA_negative_control", "DSA_positive_no_rejection", "DSA_positive_rejection")
g50084_labels <- c("DSA- / no rejection\ncontrol (n=20)", "DSA+ / no rejection\ncontrol (n=13)",
                   "DSA+ / rejection\nABMR-enriched (n=28)")
cohort_audit <- rbindlist(list(
  m50084_all[, .(profiles=.N), by=.(dataset, specimen=ifelse(compartment=="biopsy", "Kidney biopsy\n(primary)", "Whole blood\n(context only)"))],
  data.table(dataset="GSE36059", specimen="Kidney biopsy\n(validation)", profiles=381L),
  data.table(dataset="GSE275824", specimen="Kidney biopsy\n(longitudinal RCT)", profiles=63L)
))
cohort_audit[, dataset := factor(dataset, levels=c("GSE50084","GSE36059","GSE275824"))]
p_a <- ggplot(cohort_audit, aes(dataset, profiles, fill=specimen)) +
  geom_col(width=.68, color="white") +
  geom_text(aes(label=profiles), position=position_stack(vjust=.5), color="white", fontface="bold", size=3.4) +
  scale_fill_manual(values=c("Kidney biopsy\n(primary)"=unname(cols["abmr"]),
    "Whole blood\n(context only)"="#6F93A3", "Kidney biopsy\n(validation)"=unname(cols["placebo"]),
    "Kidney biopsy\n(longitudinal RCT)"=unname(cols["treated"]))) +
  labs(title="A  Cohorts and biospecimens entering the study",
       subtitle="Primary inference uses kidney-allograft biopsies; GSE50084 blood is retained but not pooled with tissue",
       x=NULL, y="Number of expression profiles", fill="Biospecimen / role") + theme_data
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
# Estimate the two nonredundant validation contrasts and display their uncertainty directly.
validation_effects <- rbindlist(lapply(selected_modules, function(mod) {
  z <- s36059[module==mod]
  rbindlist(lapply(c("non-rejecting","TCMR"), function(ref) {
    fit <- t.test(z[analysis_group=="ABMR", score], z[analysis_group==ref, score])
    data.table(module=mod, comparison=ifelse(ref=="non-rejecting",
      "Replication: ABMR vs non-rejecting", "Specificity: ABMR vs TCMR"),
      estimate=unname(fit$estimate[1]-fit$estimate[2]), lo=fit$conf.int[1],
      hi=fit$conf.int[2], p_value=fit$p.value)
  }))
}))
validation_effects[, fdr := p.adjust(p_value, method="BH"), by=comparison]
validation_effects[, module_label := factor(module_labels[module], levels=rev(module_labels[selected_modules]))]
validation_effects[, evidence := ifelse(fdr<.05,"FDR < 0.05","FDR >= 0.05")]
fwrite(validation_effects, file.path(out_dir, "GSE36059_validation_effects.csv"))
p_c <- ggplot(validation_effects, aes(estimate,module_label,color=comparison)) +
  geom_vline(xintercept=0,color="grey75",linewidth=.4) +
  geom_errorbar(aes(xmin=lo,xmax=hi),orientation="y",width=.12,position=position_dodge(.55)) +
  geom_point(aes(shape=evidence),size=2.6,position=position_dodge(.55)) +
  scale_color_manual(values=c("Replication: ABMR vs non-rejecting"=unname(cols["abmr"]),
                              "Specificity: ABMR vs TCMR"=unname(cols["tcmr"]))) +
  scale_shape_manual(values=c("FDR < 0.05"=16,"FDR >= 0.05"=1)) +
  labs(title="C  Independent cohort: replication is not the same as specificity",
       subtitle="GSE36059 kidney biopsies; right of zero means a higher module score in ABMR",
       x="ABMR minus comparator (mean difference, 95% CI)",y=NULL,
       color="Prespecified comparison",shape="Evidence") + theme_data
p_d <- make_distributions(s36059, "analysis_group", g36059_levels, g36059_labels,
                          unname(cols[c("control", "tcmr", "abmr")]),
                          "D  GSE36059: observed biopsy-level module scores",
                          "All 381 prespecified validation biopsies are displayed")

# GSE275824: visualize actual longitudinal patient measurements and randomized assignment.
m275824 <- fread(file.path(root, "data/processed/GSE275824/samples.csv"))
score_labels <- c(
  abmract = "Overall ABMR activity (ABMRact)",
  abmrifng = "IFN-gamma-inducible ABMR genes (ABMRIFNG)",
  abmrnk = "NK-expressed ABMR genes (ABMRNK)",
  abmrendo = "ABMR-associated endothelial genes (ABMREndo)")
trajectory <- melt(m275824,
                   id.vars = c("sample_accession", "patient_id", "treatment", "timepoint"),
                   measure.vars = names(score_labels), variable.name = "score_code",
                   value.name = "score")
trajectory[, published_score := factor(score_labels[as.character(score_code)],
                                       levels = unname(score_labels))]
trajectory[, timepoint := factor(timepoint, levels = c("Baseline", "Week24", "Week52"))]
trajectory[, treatment_label := factor(treatment, levels = c("Placebo", "Felzartamab"),
                                       labels = c("Placebo (control)", "Felzartamab (treated)"))]
fwrite(trajectory, file.path(out_dir, "GSE275824_patient_trajectories.csv"))
mean_trajectory <- trajectory[, .(mean = mean(score, na.rm = TRUE),
                                  se = sd(score, na.rm = TRUE) / sqrt(sum(is.finite(score)))),
                              by = .(published_score, treatment_label, timepoint)]
baseline <- trajectory[timepoint=="Baseline", .(patient_id, published_score, baseline=score)]
change <- merge(trajectory[timepoint!="Baseline"], baseline, by=c("patient_id","published_score"))
change[, change_from_baseline := score-baseline]
fwrite(change, file.path(out_dir, "GSE275824_change_from_baseline.csv"))

# GSE275824 longitudinal figure: four aligned panels using original deposited transcript scores.
# Statistical annotations are read from the executed patient-blocked mixed-effects model.
interaction_stats <- fread(file.path(root, "results/bulk/GSE275824/published_score_interactions.csv"))
panel_names <- c(
  abmract="E  Overall ABMR activity (ABMRact)",
  abmrifng="F  IFN-gamma-inducible ABMR genes (ABMRIFNG)",
  abmrnk="G  NK-expressed ABMR genes (ABMRNK)",
  abmrendo="H  ABMR-associated endothelial genes (ABMREndo)")

make_longitudinal_panel <- function(code) {
  z <- trajectory[score_code==code]
  mz <- mean_trajectory[published_score==score_labels[[code]]]
  st <- interaction_stats[score==code]
  w24 <- st[grepl("Week24",term)][1]
  w52 <- st[grepl("Week52",term)][1]
  stat_label <- sprintf(
    "Treatment x time:  W24 beta = %.2f, FDR = %.3g   |   W52 beta = %.2f, FDR = %.3g",
    w24$Value,w24$fdr,w52$Value,w52$fdr)
  ggplot(z,aes(timepoint,score,group=interaction(patient_id,treatment_label),color=treatment_label)) +
    geom_line(alpha=.24,linewidth=.48) +
    geom_point(alpha=.50,size=1.25) +
    geom_errorbar(data=mz,
      aes(x=timepoint,y=mean,ymin=mean-1.96*se,ymax=mean+1.96*se,
          color=treatment_label,group=treatment_label),
      inherit.aes=FALSE,width=.06,linewidth=.72) +
    geom_line(data=mz,aes(x=timepoint,y=mean,color=treatment_label,group=treatment_label),
      inherit.aes=FALSE,linewidth=1.15) +
    geom_point(data=mz,aes(x=timepoint,y=mean,color=treatment_label),
      inherit.aes=FALSE,size=2.5) +
    annotate("label",x=1,y=Inf,label=stat_label,hjust=0,vjust=1.15,
      size=2.7,linewidth=.18,color="grey20",fill="white") +
    scale_color_manual(values=c("Placebo (control)"=unname(cols["placebo"]),
                                "Felzartamab (treated)"=unname(cols["treated"]))) +
    scale_x_discrete(labels=c(Baseline="Baseline",Week24="Week 24",Week52="Week 52")) +
    labs(title=panel_names[[code]],x="Protocol-biopsy time point",
      y="Deposited MMDx transcript-derived score (dimensionless)",
      color="Randomized treatment arm") +
    theme_classic(base_size=9.5) +
    theme(plot.title=element_text(face="bold",size=10.5,margin=margin(b=2)),
      axis.title=element_text(face="bold",size=9.5),
      axis.text=element_text(size=8.5,color="black"),
      axis.title.x=element_text(margin=margin(t=2)),
      axis.title.y=element_text(margin=margin(r=2)),
      legend.position="bottom",
      legend.title=element_text(face="bold",size=7.5),
      legend.text=element_text(size=7),
      legend.key.width=grid::unit(10,"pt"),
      legend.background=element_rect(fill="white",color="grey45",linewidth=.3),
      legend.box.background=element_rect(fill="white",color="grey45",linewidth=.3),
      legend.margin=margin(2,3,2,3),
      plot.margin=margin(3,3,3,3))
}
longitudinal_panels <- lapply(names(panel_names),make_longitudinal_panel)
p_spaghetti <- wrap_plots(longitudinal_panels,ncol=2,guides="collect") &
  theme(legend.position="bottom",
    legend.background=element_rect(fill="white",color="grey45",linewidth=.3),
    legend.box.background=element_rect(fill="white",color="grey45",linewidth=.3))
# The longitudinal panels are embedded in Figure 1; no redundant standalone plot is exported.

# Export audit tables containing every expression-derived point used in panels B and D.
audit_scores <- rbindlist(list(s50084, s36059), fill = TRUE)
fwrite(audit_scores, file.path(out_dir, "individual_sample_module_scores.csv"))

# Figure 1 combines the four cohort panels and four longitudinal panels as A-H.
cohort_block <- (p_a | p_b) / (p_c | p_d) + plot_layout(heights=c(.86,1))
longitudinal_block <- wrap_plots(longitudinal_panels,ncol=2,guides="collect") &
  theme(legend.position="bottom")
main <- wrap_elements(full=cohort_block) / wrap_elements(full=longitudinal_block) +
  plot_layout(heights=c(1,1.55)) +
  plot_annotation(caption=paste0(
    "Panels E-H use original deposited GSE275824 values. Thin lines are patients; thick lines are arm means with approximate 95% CI. ",
    "Beta is the felzartamab-minus-placebo treatment-by-time interaction. Diebold et al., Nature Medicine 2025, doi:10.1038/s41591-025-03653-3."),
    theme=theme(plot.caption=element_text(size=7.5,color="grey25",hjust=0))) &
  theme(plot.margin=margin(2,2,2,2))

ggsave(file.path(out_dir, "dataset_results_main.png"),main,width=16,height=21,dpi=300,bg="white")
ggsave(file.path(out_dir, "dataset_results_main.pdf"),main,width=16,height=21,device="pdf")

cat("Dataset-level figures generated from", nrow(m50084), "GSE50084 biopsies,",
    nrow(m36059), "GSE36059 biopsies, and", nrow(m275824), "GSE275824 biopsies.\n")
