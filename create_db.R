#Code to generate list of lists to be fed into glm_predict.R

library(MatrixEQTL)
library(plyr)
library(BatchExperiments)
library(doParallel)
library(sqldf)
library(RSQLite)
library(reshape2)
if(Sys.info()['sysname']=="Windows"){
  snp.exploc <- "C:/Users/nknoblau/Documents/R_WS/MatrixeQTLGLM/test/snp.exp.Rdata"
  
}else{
  snp.exploc <- "/scratch/nwk2/mEQTL_ERpnc/glmEQTL/unimputed_brca/snp.exp.Rdata"
}

load(snp.exploc)
snp.exp$snps <- t(as.matrix(snp.exp$snps))
snp.exp$gene <- t(as.matrix(snp.exp$gene))
snp.exp$gene <- snp.exp$gene[rownames(snp.exp$snps),]

db <- dbConnect(drv=dbDriver("SQLite"),"2ngene_snp.db",synchronous=0)

gene.melt <- melt(snp.exp$gene)
colnames(gene.melt) <- c("Sample","Gene","Value")
dbRemoveTable(db,"genes")
dbSendQuery(db,"CREATE TABLE genes (Sample, Gene Text, Value REAL)")
dbSendQuery(db,"CREATE INDEX gs ON genes(Gene ASC,Sample ASC)")

gene.melt <- gene.melt[order(gene.melt$Gene,gene.melt$Sample),]
dbWriteTable(db,"genes",gene.melt,append=T,overwrite=F,row.names=F)

dbDisconnect(db);
