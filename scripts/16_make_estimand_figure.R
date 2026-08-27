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
theme_paper <- theme_classic(base_size=10.5) + theme(
  plot.title=element_text(face="bold", size=11.5), plot.subtitle=element_text(size=9, color="grey30"),
  strip.background=element_rect(fill="grey94", color=NA), strip.text=element_text(face="bold"),
  axis.title=element_text(face="bold"), legend.position="bottom", legend.title=element_text(face="bold"))

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
       x="Case - control module effect",y=NULL,color="Estimand",shape="Evidence") + theme_paper

# B: concordance is not sufficient; color explicitly reports the ABMR-versus-TCMR estimate.
disc <- e[contrast=="ABMR_DSApos_vs_DSApos_nonrejection",.(module,discovery=logFC,discovery_fdr=module_fdr)]
val <- e[contrast=="ABMR_vs_nonrejection",.(module,validation=logFC,validation_fdr=module_fdr)]
spec <- e[contrast=="ABMR_vs_TCMR",.(module,specificity=logFC,specificity_fdr=module_fdr)]
rep <- Reduce(function(x,y) merge(x,y,by="module"),list(disc,val,spec))[module %in% primary]
rep[, module_label:=labels[module]]
rep[, replicated:=discovery_fdr<.05 & validation_fdr<.05 & sign(discovery)==sign(validation)]
p_b <- ggplot(rep,aes(discovery,validation,color=specificity)) +
  geom_hline(yintercept=0,color="grey80") + geom_vline(xintercept=0,color="grey80") +
  geom_point(aes(shape=replicated),size=3.2) +
  geom_text(aes(label=module_label),nudge_y=.055,size=2.8,check_overlap=TRUE) +
  scale_color_gradient2(low=pal["purple"],mid="grey92",high=pal["red"],midpoint=0) +
  scale_shape_manual(values=c(`TRUE`=16,`FALSE`=1)) +
  labs(title="B  Replication is challenged for rejection specificity",
       subtitle="Color: ABMR - TCMR; purple means the program is stronger in TCMR",
       x="GSE50084: DSA+ rejection - DSA+ non-rejection",y="GSE36059: ABMR - non-rejecting",
       color="ABMR - TCMR",shape="Replicated") + theme_paper

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
       subtitle="Red: enriched in first group; blue: enriched in control/placebo direction; opacity denotes BH FDR",
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
p_d <- ggplot(d,aes(Value,score_label,color=week)) + geom_vline(xintercept=0,color="grey75") +
  geom_errorbar(aes(xmin=lo,xmax=hi),orientation="y",width=.13,position=position_dodge(.48)) +
  geom_point(aes(shape=evidence),size=2.8,position=position_dodge(.48)) +
  scale_color_manual(values=c("Week 24"=pal["teal"],"Week 52"=pal["red"])) +
  scale_shape_manual(values=c("FDR < 0.05"=16,"FDR >= 0.05"=1)) +
  labs(title="D  Randomized treatment effects are transient and program-dependent",
       subtitle="Felzartamab-minus-placebo difference-in-differences; negative values indicate relative reduction",
       x="Treatment-by-time interaction (95% CI)",y=NULL,color="Time",shape="Evidence") + theme_paper

fig <- (p_a|p_b)/(p_c|p_d) + plot_annotation(
  title="MOSAIC-Tx: qualification of a treatment-responsive rejection program",
  subtitle="Effect estimates answer prespecified biological questions; descriptive PCA is excluded from the primary result",
  caption="Associative bulk transcriptomics cannot prove antibody engagement, cell of origin, target causality, or drug efficacy.",
  theme=theme(plot.title=element_text(face="bold",size=15),plot.subtitle=element_text(size=10),plot.caption=element_text(size=8)))
ggsave(file.path(out,"estimand_results_main.png"),fig,width=16,height=12,dpi=300,bg="white")
ggsave(file.path(out,"estimand_results_main.pdf"),fig,width=16,height=12,device="pdf")
fwrite(rep,file.path(out,"cross_cohort_replication_specificity.csv"))
fwrite(g,file.path(out,"selected_hallmark_gsea_atlas.csv"))
cat("Estimand-centered biological results figure generated.\n")
