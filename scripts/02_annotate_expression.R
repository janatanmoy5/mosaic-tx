#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(jsonlite)})
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3) stop("usage: 02_annotate_expression.R <probe_matrix.csv.gz> <GPL.annot.gz> <output_dir>")
expr_path<-args[1]; annot_path<-args[2]; out<-args[3]; dir.create(out,recursive=TRUE,showWarnings=FALSE)
lines<-readLines(gzfile(annot_path),n=80); header<-grep("^ID\\t",lines)[1]
if(is.na(header)) stop("annotation header not found")
ann<-fread(cmd=paste("gzip -cd",shQuote(annot_path)),skip=header-1,select=c("ID","Gene symbol"),na.strings=c("","---"))
setnames(ann,c("probe_id","gene_symbol")); ann[,gene_symbol:=sub(" ///.*$","",gene_symbol)]
ann<-unique(ann[!is.na(gene_symbol) & gene_symbol!=""],by="probe_id")
expr<-fread(expr_path); setnames(expr,1,"probe_id"); merged<-merge(ann,expr,by="probe_id")
if(nrow(merged)<1000) stop("too few annotated probes")
sample_cols<-setdiff(names(merged),c("probe_id","gene_symbol"))
for(n in sample_cols) set(merged,j=n,value=as.numeric(merged[[n]]))
gene<-merged[,lapply(.SD,mean,na.rm=TRUE),by=gene_symbol,.SDcols=sample_cols]
fwrite(gene,file.path(out,"gene_expression.csv.gz"))
manifest<-list(schema_version="1.0",input_probes=nrow(expr),annotated_probes=nrow(merged),unique_genes=nrow(gene),aggregation="arithmetic mean of submitted normalized probe values",annotation_platform="GPL6244",model_ready=FALSE,blocking_review="clinical phenotype review")
write_json(manifest,file.path(out,"annotation_manifest.json"),pretty=TRUE,auto_unbox=TRUE)
print(manifest)
