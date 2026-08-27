#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(limma);library(fgsea);library(msigdbr);library(jsonlite)})
a<-commandArgs(trailingOnly=TRUE);expr<-fread(a[1]);s<-fread(a[2]);out<-a[3];dir.create(out,recursive=TRUE,showWarnings=FALSE)
m<-as.matrix(expr[,-1]);rownames(m)<-expr$gene_symbol;idx<-match(s$sample_accession,colnames(m));if(anyNA(idx))stop("sample mismatch");m<-m[,idx]
sets<-msigdbr(species="Homo sapiens",collection="H");paths<-split(sets$gene_symbol,sets$gs_name);summary<-list()
for(control in c("non-rejecting","TCMR")){
 keep<-s$analysis_group%in%c("ABMR",control);grp<-factor(s$analysis_group[keep],levels=c("ABMR",control));design<-model.matrix(~0+grp);colnames(design)<-make.names(levels(grp));fit<-eBayes(contrasts.fit(lmFit(m[,keep],design),makeContrasts(ABMR-get(make.names(control)),levels=design)))
 tab<-as.data.table(topTable(fit,n=Inf),keep.rownames="gene_symbol");nm<-paste0("ABMR_vs_",gsub("-","_",control));fwrite(tab,file.path(out,paste0("de_",nm,".csv")));ranks<-sort(setNames(tab$t,tab$gene_symbol),decreasing=TRUE);fg<-as.data.table(fgseaMultilevel(paths,ranks,minSize=10,maxSize=500));setorder(fg,padj);fg[,leadingEdge:=vapply(leadingEdge,paste,collapse=";",FUN.VALUE=character(1))];fwrite(fg,file.path(out,paste0("gsea_",nm,".csv")));summary[[nm]]<-list(n=sum(keep),counts=as.list(table(grp)),fdr=sum(tab$adj.P.Val<.05),top_gene=tab$gene_symbol[1],top_pathway=fg$pathway[1])
}
write_json(summary,file.path(out,"validation_manifest.json"),pretty=TRUE,auto_unbox=TRUE)
