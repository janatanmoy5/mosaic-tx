#!/usr/bin/env Rscript

# Purpose: construct the public-data ABMR signature robustness atlas used by the pilot study.
# Scientific unit: a biopsy/sample for cross-sectional cohorts and a patient for longitudinal models.
# Important: module scores are relative within each cohort and are never pooled across platforms.

# Load only declared statistical and input/output dependencies.
suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(limma)
})

# Require the repository root and output directory so execution is location-independent.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: 08_build_robustness_atlas.R <repo_root> <output_dir>")
root <- normalizePath(args[1], mustWork = TRUE)
out <- args[2]
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# Read the frozen, outcome-independent module definitions.
module_cfg <- fromJSON(file.path(root, "configs", "mechanism_modules.json"), simplifyVector = FALSE)
modules <- lapply(module_cfg$modules, function(x) unlist(x$genes, use.names = FALSE))
min_coverage <- as.numeric(module_cfg$rules$minimum_gene_coverage_fraction)

# Read a gene-by-sample matrix and collapse any duplicate symbols deterministically by their mean.
read_expression <- function(dataset) {
  path <- file.path(root, "data", "processed", dataset, "gene_expression.csv.gz")
  x <- fread(path)
  setnames(x, 1, "gene_symbol")
  x[, gene_symbol := toupper(as.character(gene_symbol))]
  sample_cols <- setdiff(names(x), "gene_symbol")
  x <- x[, lapply(.SD, mean, na.rm = TRUE), by = gene_symbol, .SDcols = sample_cols]
  mat <- as.matrix(x[, ..sample_cols])
  rownames(mat) <- x$gene_symbol
  storage.mode(mat) <- "double"
  mat
}

# Convert each module into the mean of gene-wise z scores calculated within one cohort.
# This avoids merging measurement scales from three different microarray platforms.
score_modules <- function(mat, dataset) {
  score_rows <- list()
  coverage_rows <- list()
  for (module_name in names(modules)) {
    declared <- unique(toupper(modules[[module_name]]))
    present <- intersect(declared, rownames(mat))
    coverage <- length(present) / length(declared)
    coverage_rows[[module_name]] <- data.table(
      dataset = dataset,
      module = module_name,
      declared_genes = length(declared),
      observed_genes = length(present),
      coverage_fraction = coverage,
      observed_gene_symbols = paste(present, collapse = ";")
    )
    if (coverage < min_coverage) next
    block <- mat[present, , drop = FALSE]
    gene_sd <- apply(block, 1, sd, na.rm = TRUE)
    block <- block[is.finite(gene_sd) & gene_sd > 0, , drop = FALSE]
    if (!nrow(block)) next
    z <- t(scale(t(block)))
    score_rows[[module_name]] <- colMeans(z, na.rm = TRUE)
  }
  if (!length(score_rows)) stop("no module passed coverage for ", dataset)
  list(
    scores = do.call(rbind, score_rows),
    coverage = rbindlist(coverage_rows, use.names = TRUE)
  )
}

# Fit a two-group limma contrast and return one row per prespecified module.
fit_pair <- function(scores, metadata, group_column, case, control, dataset, contrast_name) {
  keep <- metadata[[group_column]] %in% c(case, control)
  d <- metadata[keep]
  idx <- match(d$sample_accession, colnames(scores))
  if (anyNA(idx)) stop("sample/matrix mismatch in ", dataset, " ", contrast_name)
  y <- scores[, idx, drop = FALSE]
  group <- factor(d[[group_column]], levels = c(control, case))
  design <- model.matrix(~ group)
  fit <- eBayes(lmFit(y, design))
  tab <- as.data.table(topTable(fit, coef = 2, number = Inf, sort.by = "none"), keep.rownames = "module")
  tab[, `:=`(
    dataset = dataset,
    contrast = contrast_name,
    case_group = case,
    control_group = control,
    n_case = sum(group == case),
    n_control = sum(group == control)
  )]
  tab
}

# Fit placebo-adjusted longitudinal module changes while retaining within-patient correlation.
fit_treatment <- function(scores, metadata) {
  idx <- match(metadata$sample_accession, colnames(scores))
  if (anyNA(idx)) stop("sample/matrix mismatch in GSE275824")
  y <- scores[, idx, drop = FALSE]
  group <- factor(paste(metadata$timepoint, metadata$treatment, sep = "_"))
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  block <- factor(metadata$patient_id)
  corfit <- duplicateCorrelation(y, design, block = block)
  fit <- lmFit(y, design, block = block, correlation = corfit$consensus)
  contrast_matrix <- makeContrasts(
    W24 = (Week24_Felzartamab - Baseline_Felzartamab) - (Week24_Placebo - Baseline_Placebo),
    W52 = (Week52_Felzartamab - Baseline_Felzartamab) - (Week52_Placebo - Baseline_Placebo),
    levels = design
  )
  fit <- eBayes(contrasts.fit(fit, contrast_matrix))
  rows <- lapply(seq_len(ncol(contrast_matrix)), function(i) {
    tab <- as.data.table(topTable(fit, coef = i, number = Inf, sort.by = "none"), keep.rownames = "module")
    tab[, `:=`(
      dataset = "GSE275824",
      contrast = paste0("treatment_difference_in_differences_", colnames(contrast_matrix)[i]),
      case_group = "Felzartamab",
      control_group = "Placebo",
      n_case = uniqueN(metadata[treatment == "Felzartamab"]$patient_id),
      n_control = uniqueN(metadata[treatment == "Placebo"]$patient_id)
    )]
    tab
  })
  list(table = rbindlist(rows), consensus_correlation = unname(corfit$consensus))
}

# Compile all three public cohorts independently.
datasets <- c("GSE50084", "GSE36059", "GSE275824")
compiled <- lapply(datasets, function(dataset) {
  mat <- read_expression(dataset)
  scored <- score_modules(mat, dataset)
  meta <- fread(file.path(root, "data", "processed", dataset, "samples.csv"))
  list(dataset = dataset, scores = scored$scores, coverage = scored$coverage, metadata = meta)
})
names(compiled) <- datasets

# Declare cross-sectional contrasts, including the difficult DSA-positive/no-rejection reference.
g50084_meta <- compiled$GSE50084$metadata[compartment == "biopsy"]
effects <- list(
  fit_pair(compiled$GSE50084$scores, g50084_meta, "analysis_group",
           "DSA_positive_rejection", "DSA_negative_control", "GSE50084", "ABMR_DSApos_vs_DSAneg_nonrejection"),
  fit_pair(compiled$GSE50084$scores, g50084_meta, "analysis_group",
           "DSA_positive_rejection", "DSA_positive_no_rejection", "GSE50084", "ABMR_DSApos_vs_DSApos_nonrejection"),
  fit_pair(compiled$GSE50084$scores, g50084_meta, "analysis_group",
           "DSA_positive_no_rejection", "DSA_negative_control", "GSE50084", "DSApos_nonrejection_vs_DSAneg_nonrejection"),
  fit_pair(compiled$GSE36059$scores, compiled$GSE36059$metadata, "analysis_group",
           "ABMR", "non-rejecting", "GSE36059", "ABMR_vs_nonrejection"),
  fit_pair(compiled$GSE36059$scores, compiled$GSE36059$metadata, "analysis_group",
           "ABMR", "TCMR", "GSE36059", "ABMR_vs_TCMR")
)

# Add the randomized longitudinal contrasts and control multiplicity within each contrast family.
treatment <- fit_treatment(compiled$GSE275824$scores, compiled$GSE275824$metadata)
effects <- rbindlist(c(effects, list(treatment$table)), fill = TRUE)
effects[, module_fdr := p.adjust(P.Value, method = "BH"), by = .(dataset, contrast)]
setcolorder(effects, c("dataset", "contrast", "module", "case_group", "control_group",
                       "n_case", "n_control", "logFC", "t", "P.Value", "module_fdr"))
setorder(effects, contrast, module)

# Translate the effect table into a conservative, transparent evidence grade.
wide <- dcast(effects, module ~ contrast, value.var = c("logFC", "module_fdr"))
discovery <- "ABMR_DSApos_vs_DSAneg_nonrejection"
validation <- "ABMR_vs_nonrejection"
specificity <- "ABMR_vs_TCMR"
wide[, direction_concordant := sign(get(paste0("logFC_", discovery))) == sign(get(paste0("logFC_", validation)))]
wide[, replicated_fdr := get(paste0("module_fdr_", discovery)) < 0.05 & get(paste0("module_fdr_", validation)) < 0.05]
wide[, tcmr_specific_fdr := get(paste0("module_fdr_", specificity)) < 0.05]
wide[, evidence_tier := fifelse(direction_concordant & replicated_fdr & tcmr_specific_fdr, "A",
                         fifelse(direction_concordant & replicated_fdr, "B",
                           fifelse(direction_concordant &
                             (get(paste0("module_fdr_", discovery)) < 0.05 |
                              get(paste0("module_fdr_", validation)) < 0.05), "C", "D")))]
wide[, evidence_definition := "A=replicated+TCMR-specific; B=replicated not TCMR-specific; C=partial; D=unsupported"]

# Save auditable intermediate and final outputs rather than only a rendered figure.
coverage <- rbindlist(lapply(compiled, `[[`, "coverage"))
fwrite(coverage, file.path(out, "module_coverage.csv"))
fwrite(effects, file.path(out, "module_effects.csv"))
fwrite(wide, file.path(out, "robustness_atlas.csv"))

# Write a machine-readable run manifest with methods, sample sizes, and claim boundaries.
manifest <- list(
  schema_version = "1.0",
  generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  datasets = datasets,
  modules = names(modules),
  contrasts = unique(effects$contrast),
  treatment_consensus_correlation = treatment$consensus_correlation,
  methods = c("within-cohort gene z-score module means", "limma moderated contrasts",
              "patient-blocked duplicateCorrelation longitudinal difference-in-differences", "BH FDR per contrast"),
  limitations = c("retrospective public cohorts", "module effects are association not causation",
                  "bulk modules cannot establish cell lineage", "GSE275824 contains 22 randomized patients",
                  "no clinical model or drug efficacy claim")
)
write_json(manifest, file.path(out, "run_manifest.json"), pretty = TRUE, auto_unbox = TRUE)
print(manifest)
