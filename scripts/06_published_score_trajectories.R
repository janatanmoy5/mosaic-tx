#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(nlme);library(jsonlite)})
a<-commandArgs(trailingOnly=TRUE);if(length(a)!=2)stop("usage: samples.csv output_dir")
s<-fread(a[1]);out<-a[2];dir.create(out,recursive=TRUE,showWarnings=FALSE)
s[,patient_id:=factor(patient_id)];s[,treatment:=factor(treatment,levels=c("Placebo","Felzartamab"))];s[,timepoint:=factor(timepoint,levels=c("Baseline","Week24","Week52"))]
scores<-c("abmract","abmrifng","abmrnk","abmrendo","nkb","endat");rows<-list()
for(score in scores){
 d<-s[,.(patient_id,treatment,timepoint,value=as.numeric(get(score)))][is.finite(value)]
 fit<-lme(value~treatment*timepoint,random=~1|patient_id,data=d,method="REML",na.action=na.omit)
 tab<-as.data.table(summary(fit)$tTable,keep.rownames="term");tab[,score:=score];rows[[score]]<-tab
}
res<-rbindlist(rows,use.names=TRUE,fill=TRUE);setcolorder(res,c("score","term",setdiff(names(res),c("score","term"))));fwrite(res,file.path(out,"published_score_mixed_models.csv"))
interactions<-res[grepl("treatmentFelzartamab:timepoint",term)];interactions[,fdr:=p.adjust(`p-value`,method="BH")];fwrite(interactions,file.path(out,"published_score_interactions.csv"))
manifest<-list(method="linear mixed effects; patient random intercept; treatment-by-time interaction; BH across declared scores/timepoints",n_patients=uniqueN(s$patient_id),n_biopsies=nrow(s),scores=scores,significant_interactions_fdr_0_05=sum(interactions$fdr<.05));write_json(manifest,file.path(out,"published_score_manifest.json"),pretty=TRUE,auto_unbox=TRUE);print(manifest)
