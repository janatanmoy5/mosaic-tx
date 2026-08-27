#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(AnnotationDbi);library(jsonlite)})
a<-commandArgs(trailingOnly=TRUE);if(length(a)!=4)stop("usage: matrix samples platform out")
expr<-fread(a[1]);setnames(expr,1,"probe_id"); platform<-a[3]; out<-a[4];dir.create(out,recursive=TRUE,showWarnings=FALSE)
pkg<-switch(platform,GPL6244="hugene10sttranscriptcluster.db",GPL570="hgu133plus2.db",GPL13667="hgu219.db",stop("unsupported platform"))
suppressPackageStartupMessages(library(pkg,character.only=TRUE));db<-get(pkg)
map<-AnnotationDbi::mapIds(db,keys=as.character(expr$probe_id),column="SYMBOL",keytype="PROBEID",multiVals="first")
expr[,gene_symbol:=unname(map[as.character(probe_id)])]; sample_cols<-setdiff(names(expr),c("probe_id","gene_symbol")); annotated<-expr[!is.na(gene_symbol)&gene_symbol!=""]
for(n in sample_cols)set(annotated,j=n,value=as.numeric(annotated[[n]]))
gene<-annotated[,lapply(.SD,mean,na.rm=TRUE),by=gene_symbol,.SDcols=sample_cols];fwrite(gene,file.path(out,"gene_expression.csv.gz"))
manifest<-list(platform=platform,annotation_package=pkg,input_probes=nrow(expr),annotated_probes=nrow(annotated),unique_genes=nrow(gene),aggregation="mean",model_ready=TRUE);write_json(manifest,file.path(out,"annotation_manifest.json"),pretty=TRUE,auto_unbox=TRUE);print(manifest)
