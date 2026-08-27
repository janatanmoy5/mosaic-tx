#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(limma);library(fgsea);library(msigdbr);library(jsonlite)})
a<-commandArgs(trailingOnly=TRUE);x<-fread(a[1]);s<-fread(a[2]);out<-a[3];dir.create(out,recursive=TRUE,showWarnings=FALSE)
m<-as.matrix(x[,-1]);rownames(m)<-x$gene_symbol;idx<-match(s$sample_accession,colnames(m));if(anyNA(idx))stop("sample mismatch");m<-m[,idx]
group<-factor(paste(s$timepoint,s$treatment,sep="_"));design<-model.matrix(~0+group);colnames(design)<-levels(group);block<-factor(s$patient_id)
corfit<-duplicateCorrelation(m,design,block=block);fit<-lmFit(m,design,block=block,correlation=corfit$consensus)
cons<-makeContrasts(W24=(Week24_Felzartamab-Baseline_Felzartamab)-(Week24_Placebo-Baseline_Placebo),W52=(Week52_Felzartamab-Baseline_Felzartamab)-(Week52_Placebo-Baseline_Placebo),levels=design);fit<-eBayes(contrasts.fit(fit,cons))
sets<-msigdbr(species="Homo sapiens",collection="H");paths<-split(sets$gene_symbol,sets$gs_name);summary<-list(consensus_correlation=corfit$consensus)
for(i in 1:2){nm<-colnames(cons)[i];tab<-as.data.table(topTable(fit,coef=i,n=Inf),keep.rownames="gene_symbol");fwrite(tab,file.path(out,paste0("de_difference_in_differences_",nm,".csv")));ranks<-sort(setNames(tab$t,tab$gene_symbol),decreasing=TRUE);fg<-as.data.table(fgseaMultilevel(paths,ranks,minSize=10,maxSize=500));setorder(fg,padj);fg[,leadingEdge:=vapply(leadingEdge,paste,collapse=";",FUN.VALUE=character(1))];fwrite(fg,file.path(out,paste0("gsea_difference_in_differences_",nm,".csv")));summary[[nm]]<-list(fdr=sum(tab$adj.P.Val<.05),top_gene=tab$gene_symbol[1],top_pathway=fg$pathway[1],pathways_fdr=sum(fg$padj<.05,na.rm=TRUE))}
summary$design<-"patient-blocked limma; (follow-up - baseline) felzartamab minus placebo";summary$n_patients<-length(unique(s$patient_id));summary$n_biopsies<-nrow(s);write_json(summary,file.path(out,"treatment_manifest.json"),pretty=TRUE,auto_unbox=TRUE)
