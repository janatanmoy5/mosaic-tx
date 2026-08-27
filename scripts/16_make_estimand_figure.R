#!/usr/bin/env Rscript

# Build the primary biological-results figure from saved, executed statistical estimates.
# Unlike PCA, every panel corresponds to a prespecified biological estimand.

suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
a <- commandArgs(trailingOnly = TRUE)
if (length(a) != 3) stop("usage: 16_make_estimand_figure.R <repo_root> <atlas_dir> <output_dir>")
root <- normalizePath(a[1], mustWork = TRUE); atlas <- a[2]; out <- a[3]
dir.create(out, recursive = TRUE, showWarnings = FALSE)

pal <- c(red="#B64B4B", gold="#D39C3F", navy="#24476B", purple="#6F5AA8",
         teal="#16817A", grey="#707B83")
theme_paper <- theme_classic(base_size=9.5) + theme(
  plot.title=element_text(face="bold",size=10.5,margin=margin(b=2)),
  plot.subtitle=element_text(size=8,color="grey30",margin=margin(b=2)),
  strip.background=element_rect(fill="grey94",color=NA),
  strip.text=element_text(face="bold",size=8.5),
  axis.title=element_text(face="bold",size=9.5),
  axis.title.x=element_text(margin=margin(t=2)),
  axis.title.y=element_text(margin=margin(r=0)),
  axis.text=element_text(size=8.5,color="black"),
  axis.text.y=element_text(margin=margin(r=0)),
  legend.position="bottom",
  legend.title=element_text(face="bold",size=7.5),
  legend.text=element_text(size=7),
  legend.key.height=grid::unit(8,"pt"),
  legend.key.width=grid::unit(9,"pt"),
  legend.spacing.x=grid::unit(2,"pt"),
  legend.background=element_rect(fill="white",color="grey45",linewidth=.3),
  legend.box.background=element_rect(fill="white",color="grey45",linewidth=.3),
  legend.margin=margin(2,3,2,3),
  plot.margin=margin(2,2,2,2))

labels <- c(nk_identity="NK identity", fcgr3a_adcc_effector="FCGR3A/Fc signaling",
  cytotoxicity="Cytotoxicity", ifng_chemokine="IFN-gamma/chemokines",
  endothelial_identity="Endothelial identity", endothelial_activation_injury="Endothelial activation",
  tissue_injury_repair="Tissue injury/repair", negative_control_t_cell="T-cell control",
  negative_control_monocyte="Monocyte control")
primary <- c("nk_identity","fcgr3a_adcc_effector","cytotoxicity","ifng_chemokine",
             "endothelial_activation_injury","tissue_injury_repair")
e <- fread(file.path(atlas,"module_effects.csv"))
e[, se := fifelse(is.finite(t) & abs(t)>0, abs(logFC/t), NA_real_)]
e[, `:=`(lo=logFC-1.96*se, hi=logFC+1.96*se, module_label=labels[module],
         evidence=ifelse(module_fdr<0.05,"FDR < 0.05","FDR >= 0.05"))]

# A: the three contrasts required to distinguish DSA exposure from damaging rejection.
a_map <- c(DSApos_nonrejection_vs_DSAneg_nonrejection="DSA exposure without rejection",
  ABMR_DSApos_vs_DSApos_nonrejection="Rejection among DSA+ grafts",
  ABMR_DSApos_vs_DSAneg_nonrejection="DSA+ rejection vs DSA- control")
pa <- e[dataset=="GSE50084" & module %in% primary & contrast %in% names(a_map)]
pa[, comparison := factor(a_map[contrast], levels=unname(a_map))]
pa[, module_label := factor(module_label, levels=rev(labels[primary]))]
p_a <- ggplot(pa, aes(logFC,module_label,color=comparison)) + geom_vline(xintercept=0,color="grey75") +
  geom_errorbar(aes(xmin=lo,xmax=hi),orientation="y",width=.12,position=position_dodge(.58)) +
  geom_point(aes(shape=evidence),size=2.5,position=position_dodge(.58)) +
  scale_color_manual(values=unname(pal[c("gold","red","navy")])) +
  scale_shape_manual(values=c("FDR < 0.05"=16,"FDR >= 0.05"=1)) +
  labs(title="A  DSA exposure is separated from damaging rejection",
       subtitle="GSE50084, 61 primary biopsies; standardized module effects with approximate 95% CI",
       x="Case - control module effect",y=NULL,color="Estimand",shape="Evidence") +
  guides(color=guide_legend(nrow=2,byrow=TRUE),shape=guide_legend(nrow=1)) + theme_paper

# B: show the three clinically distinct comparisons directly; do not encode specificity as color.
disc <- e[contrast=="ABMR_DSApos_vs_DSApos_nonrejection",.(module,estimate=logFC,lo,hi,module_fdr)]
disc[, comparison := "Discovery: DSA+ rejection vs DSA+ no rejection"]
val <- e[contrast=="ABMR_vs_nonrejection",.(module,estimate=logFC,lo,hi,module_fdr)]
val[, comparison := "Validation: ABMR vs non-rejecting"]
spec <- e[contrast=="ABMR_vs_TCMR",.(module,estimate=logFC,lo,hi,module_fdr)]
spec[, comparison := "Specificity: ABMR vs TCMR"]
rep_long <- rbindlist(list(disc,val,spec))[module %in% primary]
rep_long[, module_label := factor(labels[module], levels=rev(labels[primary]))]
rep_long[, comparison := factor(comparison, levels=c(
  "Discovery: DSA+ rejection vs DSA+ no rejection",
  "Validation: ABMR vs non-rejecting",
  "Specificity: ABMR vs TCMR"))]
rep_long[, evidence := ifelse(module_fdr<.05,"FDR < 0.05","FDR >= 0.05")]
p_b <- ggplot(rep_long,aes(estimate,module_label,color=comparison)) +
  geom_vline(xintercept=0,color="grey75") +
  geom_errorbar(aes(xmin=lo,xmax=hi),orientation="y",width=.12,position=position_dodge(.64)) +
  geom_point(aes(shape=evidence),size=2.5,position=position_dodge(.64)) +
  scale_color_manual(values=c(unname(pal["teal"]),unname(pal["navy"]),unname(pal["purple"]))) +
  scale_shape_manual(values=c("FDR < 0.05"=16,"FDR >= 0.05"=1)) +
  labs(title="B  Replication and ABMR specificity are separate tests",
       subtitle="Right of zero: higher in the first named group; ABMR-vs-TCMR is the required disease-control comparison",
       x="Standardized module difference (95% CI)",y=NULL,color="Comparison",shape="Evidence") + theme_paper
rep <- dcast(rep_long, module + module_label ~ comparison, value.var="estimate")

# C: broader pathway context from all executed Hallmark GSEA contrasts.
gsea_files <- list(
  "G50084: rejection vs DSA- NR"="results/bulk/GSE50084/gsea_hallmark_rejection_vs_dsa_negative_control.csv",
  "G50084: rejection vs DSA+ NR"="results/bulk/GSE50084/gsea_hallmark_rejection_vs_dsa_positive_control.csv",
  "G36059: ABMR vs NR"="results/bulk/GSE36059/gsea_ABMR_vs_non_rejecting.csv",
  "G36059: ABMR vs TCMR"="results/bulk/GSE36059/gsea_ABMR_vs_TCMR.csv",
  "G275824: treatment DID W24"="results/bulk/GSE275824/gsea_difference_in_differences_W24.csv",
  "G275824: treatment DID W52"="results/bulk/GSE275824/gsea_difference_in_differences_W52.csv")
g <- rbindlist(lapply(names(gsea_files),function(nm){x<-fread(file.path(root,gsea_files[[nm]]));x[,contrast:=nm];x}),fill=TRUE)
paths <- c("HALLMARK_ALLOGRAFT_REJECTION","HALLMARK_INTERFERON_GAMMA_RESPONSE",
           "HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_TNFA_SIGNALING_VIA_NFKB",
           "HALLMARK_COMPLEMENT","HALLMARK_ANGIOGENESIS","HALLMARK_HYPOXIA",
           "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
g <- g[pathway %in% paths]
g[, pathway_label:=factor(gsub("_"," ",sub("HALLMARK_","",pathway)),levels=rev(gsub("_"," ",sub("HALLMARK_","",paths))))]
g[, contrast:=factor(contrast,levels=names(gsea_files))]
p_c <- ggplot(g,aes(contrast,pathway_label,fill=NES,alpha=padj<.05)) + geom_tile(color="white",linewidth=.45) +
  scale_fill_gradient2(low=pal["navy"],mid="white",high=pal["red"],midpoint=0) +
  scale_alpha_manual(values=c(`TRUE`=1,`FALSE`=.32)) +
  labs(title="C  Pathway direction across discovery, specificity and treatment",
       subtitle="Each column names the numerator and comparator; red favors the first group, blue favors the second",
       x=NULL,y=NULL,fill="Hallmark NES",alpha="FDR < 0.05") + theme_paper +
  theme(axis.text.x=element_text(angle=32,hjust=1,size=8))

# D: randomized interaction estimates, not raw before/after differences.
d <- fread(file.path(root,"results/bulk/GSE275824/published_score_interactions.csv"))
d <- d[score %in% c("abmract","abmrifng","abmrnk","abmrendo")]
d[, `:=`(week=factor(ifelse(grepl("Week24",term),"Week 24","Week 52"),levels=c("Week 24","Week 52")),
          score_label=factor(c(abmract="ABMR activity",abmrifng="IFN-gamma",abmrnk="NK",abmrendo="Endothelium")[score],
                             levels=rev(c("ABMR activity","IFN-gamma","NK","Endothelium"))),
          lo=Value-1.96*Std.Error,hi=Value+1.96*Std.Error,
          evidence=ifelse(fdr<.05,"FDR < 0.05","FDR >= 0.05"))]
# E: platform-specific unsupervised PCA provides descriptive sample structure, not inference.
make_pca_panel <- function(dataset, keep_groups, group_labels, group_colors, panel_title) {
  meta <- fread(file.path(root,"data/processed",dataset,"samples.csv"))[
    analysis_group %in% names(group_labels)]
  expr <- fread(file.path(root,"data/processed",dataset,"gene_expression.csv.gz"))
  setnames(expr,1L,"gene_symbol")
  sample_ids <- intersect(meta$sample_accession,names(expr))
  meta <- meta[match(sample_ids,sample_accession)]
  mat <- as.matrix(expr[,..sample_ids]); storage.mode(mat) <- "double"
  gene_var <- apply(mat,1L,var,na.rm=TRUE)
  keep <- head(order(gene_var,decreasing=TRUE,na.last=NA),min(1000L,sum(is.finite(gene_var))))
  fit <- prcomp(t(mat[keep,,drop=FALSE]),center=TRUE,scale.=TRUE)
  variance <- 100*fit$sdev^2/sum(fit$sdev^2)
  pc <- data.table(sample_accession=rownames(fit$x),PC1=fit$x[,1],PC2=fit$x[,2])
  pc <- meta[pc,on="sample_accession"]
  pc[,group_label:=factor(group_labels[analysis_group],levels=unname(group_labels))]
  centroids <- pc[,.(PC1=mean(PC1),PC2=mean(PC2)),by=group_label]
  fwrite(pc,file.path(out,paste0(dataset,"_descriptive_pca.csv")))
  ggplot(pc,aes(PC1,PC2,color=group_label)) +
    geom_hline(yintercept=0,color="grey90",linewidth=.3) +
    geom_vline(xintercept=0,color="grey90",linewidth=.3) +
    stat_ellipse(aes(group=group_label),level=.80,linetype=2,linewidth=.55,show.legend=FALSE) +
    geom_point(size=1.7,alpha=.72) +
    geom_point(data=centroids,aes(PC1,PC2,fill=group_label),shape=21,size=3.2,
               color="white",stroke=.7,inherit.aes=FALSE) +
    scale_color_manual(values=setNames(group_colors,unname(group_labels)),drop=FALSE) +
    scale_fill_manual(values=setNames(group_colors,unname(group_labels)),drop=FALSE,guide="none") +
    labs(title=panel_title,
         subtitle="Unsupervised PCA of 1,000 most-variable genes; dashed contours enclose 80% of each labeled group",
         x=sprintf("PC1 (%.1f%% variance)",variance[1]),
         y=sprintf("PC2 (%.1f%% variance)",variance[2]),color="Clinical group") +
    guides(color=guide_legend(nrow=1,byrow=TRUE)) +
    theme_paper + theme(legend.position="bottom",
      axis.title.y=element_text(margin=margin(r=0)),
      plot.margin=margin(2,2,2,2))
}
p_e1 <- make_pca_panel("GSE50084",
  c("DSA_negative_control","DSA_positive_no_rejection","DSA_positive_rejection"),
  c(DSA_negative_control="DSA- / no rejection",
    DSA_positive_no_rejection="DSA+ / no rejection",
    DSA_positive_rejection="DSA+ / rejection"),
  c("#707B83","#D39C3F","#B64B4B"),
  "E  Discovery-cohort structure: GSE50084 kidney biopsies")
p_e2 <- make_pca_panel("GSE36059",c("non-rejecting","TCMR","ABMR"),
  c("non-rejecting"="Non-rejecting",TCMR="TCMR",ABMR="ABMR"),
  c("#707B83","#6F5AA8","#B64B4B"),
  "F  Validation-cohort structure: GSE36059 kidney biopsies")

p_d <- ggplot(d,aes(Value,score_label,color=week)) + geom_vline(xintercept=0,color="grey75") +
  geom_errorbar(aes(xmin=lo,xmax=hi),orientation="y",width=.13,position=position_dodge(.48)) +
  geom_point(aes(shape=evidence),size=2.8,position=position_dodge(.48)) +
  scale_color_manual(values=c("Week 24"=unname(pal["teal"]),"Week 52"=unname(pal["red"]))) +
  scale_shape_manual(values=c("FDR < 0.05"=16,"FDR >= 0.05"=1)) +
  labs(title="D  Randomized treatment effects are transient and program-dependent",
       subtitle="Felzartamab-minus-placebo difference-in-differences; negative values indicate relative reduction",
       x="Treatment-by-time interaction (95% CI)",y=NULL,color="Time",shape="Evidence") + theme_paper

pca_row <- wrap_elements(full=p_e1) | wrap_elements(full=p_e2)
fig <- (p_a|p_b)/(p_c|p_d)/pca_row +
  plot_layout(heights=c(1,1,.92),widths=c(1,1)) + plot_annotation(
  caption="Associative bulk transcriptomics cannot prove antibody engagement, cell of origin, target causality, or drug efficacy.",
  theme=theme(plot.caption=element_text(size=8,hjust=0)))
ggsave(file.path(out,"estimand_results_main.png"),fig,width=16,height=17,dpi=300,bg="white")
ggsave(file.path(out,"estimand_results_main.pdf"),fig,width=16,height=17,device="pdf")
fwrite(rep,file.path(out,"cross_cohort_replication_specificity.csv"))
fwrite(g,file.path(out,"selected_hallmark_gsea_atlas.csv"))
cat("Estimand-centered biological results figure generated.\n")
