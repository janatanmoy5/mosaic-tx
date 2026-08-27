#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(limma);library(fgsea);library(msigdbr);library(ggplot2);library(jsonlite)})
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=3) stop("usage: 03_bulk_de_gsea.R <gene_expression.csv.gz> <samples.csv> <output_dir>")
expr<-fread(args[1]); samples<-fread(args[2]); out<-args[3]; dir.create(out,recursive=TRUE,showWarnings=FALSE)
mat<-as.matrix(expr[,-1]); rownames(mat)<-expr$gene_symbol
samples<-samples[compartment=="biopsy"]; idx<-match(samples$sample_accession,colnames(mat)); if(anyNA(idx)) stop("sample/matrix mismatch"); mat<-mat[,idx,drop=FALSE]
contrasts<-list(rejection_vs_dsa_positive_control=c("DSA_positive_rejection","DSA_positive_no_rejection"),rejection_vs_dsa_negative_control=c("DSA_positive_rejection","DSA_negative_control"))
sets<-msigdbr(species="Homo sapiens",collection="H"); pathways<-split(sets$gene_symbol,sets$gs_name)
summary<-list()
for(nm in names(contrasts)){
 pair<-contrasts[[nm]]; keep<-samples$analysis_group %in% pair; y<-mat[,keep,drop=FALSE]; group<-factor(samples$analysis_group[keep],levels=pair)
 design<-model.matrix(~0+group); colnames(design)<-make.names(levels(group)); fit<-eBayes(contrasts.fit(lmFit(y,design),makeContrasts(contrasts=paste0(make.names(pair[1]),"-",make.names(pair[2])),levels=design)))
 tab<-as.data.table(topTable(fit,number=Inf,sort.by="P"),keep.rownames="gene_symbol"); fwrite(tab,file.path(out,paste0("de_",nm,".csv")))
 ranks<-tab$t; names(ranks)<-tab$gene_symbol; ranks<-sort(ranks[is.finite(ranks)],decreasing=TRUE)
 fg<-as.data.table(fgseaMultilevel(pathways,ranks,minSize=10,maxSize=500)); setorder(fg,padj); fg[,leadingEdge:=vapply(leadingEdge,paste,collapse=";",FUN.VALUE=character(1))]; fwrite(fg,file.path(out,paste0("gsea_hallmark_",nm,".csv")))
 tab[,significant:=adj.P.Val<0.05 & abs(logFC)>=0.5]
 v<-ggplot(tab,aes(logFC,-log10(pmax(adj.P.Val,1e-300)),color=significant))+geom_point(alpha=.55,size=1)+scale_color_manual(values=c("grey65","#087E8B"))+theme_minimal()+labs(title=nm,x="log2 fold change",y="-log10 FDR")
 ggsave(file.path(out,paste0("volcano_",nm,".png")),v,width=7,height=5,dpi=180)
 summary[[nm]]<-list(n_samples=ncol(y),group_counts=as.list(table(group)),fdr_0_05=sum(tab$adj.P.Val<.05),fdr_and_abs_logfc_0_5=sum(tab$significant),top_gene=tab$gene_symbol[1],top_hallmark=if(nrow(fg)) fg$pathway[1] else NA)
}
write_json(list(schema_version="1.0",method="limma moderated linear model; BH FDR; fgseaMultilevel Hallmark",contrasts=summary,limitations=c("historical submitter labels","no sample-level covariate adjustment available in GSE50084 metadata","discovery analysis only")),file.path(out,"analysis_manifest.json"),pretty=TRUE,auto_unbox=TRUE)
