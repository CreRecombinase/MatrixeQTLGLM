#Code to generate list of lists to be fed into glm_predict.R

library(MatrixEQTL)
library(plyr)
library(doParallel)
library(sqldf)
library(RSQLite)
library(RSQLite.extfuns)
library(BatchExperiments)
library(reshape2)
if(Sys.info()['sysname']=="Windows"){
  snp.exploc <- "C:/Users/nknoblau/Documents/R_WS/MatrixeQTLGLM/test/snp.exp.Rdata"
  registerDoParallel(5)
  eqtl.base <- "C:/Users/nknoblau/Documents/R_WS/MatrixeQTLGLM/test2/unimputed_brca_trans"
  snp.melt.file <- "D:/snp_melt.txt"
  
  
  
}else{
  snp.exploc <- "/scratch/nwk2/mEQTL_ERpnc/glmEQTL/unimputed_brca/snp.exp.Rdata"
  snp.melt.file <- "/scratch/nwk2/mEQTL_ERpnc/glmEQTL/unimputed_brca/snp_melt.txt"
  
}
eqtl.files <- paste0(eqtl.base,1:57,".txt")
sample.num <- 513
kfold <- 57



train.indices <- chunk(rep(1:sample.num,kfold),n.chunks=kfold)
test.indices <- chunk(1:sample.num,chunk.size=sample.num/kfold)
train.indices <- mapply(FUN=function(x,y)x[-y],train.indices,test.indices,SIMPLIFY=F)




aeqtls <- ldply(eqtl.files,function(x){
  tqtl <- read.csv.sql(x,header=T,eol="\n",sep="\t")
  knum <- gsub(".+trans([0-9]+).txt","\\1",x)
  tqtl$knum <- knum
  return(tqtl)
}
)

fsnp.genes <- split(aeqtls,aeqtls$knum)
snp.genes <- unlist(lapply(fsnp.genes,function(x)split(x$SNP,x$gene)),recursive=F)
rm(aeqtls,fsnp.genes)

load(snp.exploc)
snp.exp$snps <- t(as.matrix(snp.exp$snps))
snp.exp$gene <- t(as.matrix(snp.exp$gene))
snp.exp$gene <- snp.exp$gene[rownames(snp.exp$snps),]
gc()

db <- dbConnect(dbDriver("SQLite"),dbname="D:/gene_snp.db",loadable.extensions=T)


test.df <- do.call("rbind",lapply(1:length(test.indices),function(x){
  data.frame(Sample=rownames(snp.exp$gene)[test.indices[[x]]],Kfold=x,Index=test.indices[[x]])
}))

train.df <-test.df <- do.call("rbind",lapply(1:length(train.indices),function(x){
  data.frame(Sample=rownames(snp.exp$gene)[train.indices[[x]]],Kfold=x,Index=train.indices[[x]])
}))

dbWriteTable(db,name="testSamples",test.df,row.names=F,overwrite=T,append=F)
dbWriteTable(db,name="trainSamples",train.df,row.names=F,overwrite=T,append=F)


gene.melt <- melt(snp.exp$gene,id.vars=rownames(snp.exp$gene),variable.name="gene")
colnames(gene.melt) <- c("Sample","Gene","Value")
dbWriteTable(db,"gene",gene.melt,row.names=F,overwrite=T,append=F)
rm(gene.melt)

snpcols <- length(colnames(snp.exp$snps))
snpcol.chunks <- chunk(seq(snpcols),chunk.size=20000)


for(i in 1:length(snpcol.chunks)){
  tsnp <- melt(snp.exp$snps[,snpcol.chunks[[i]]])
  colnames(tsnp)<- c("Sample","Snps","Value")
  write.table(tsnp,file=snp.melt.file,append=T,quote=F,sep="\t",row.names=F,col.names=ifelse(i==1,T,F))
  #dbWriteTable(db,name="Snps",tsnp,row.names=F,append=T,overwrite=F)
}


dbSendQuery(conn=db,statement="CREATE INDEX gs ON gene(Gene,Sample)")
dbSendQuery(conn=db,statement="CREATE INDEX ss ON Snps(Snps,Sample)")



system.time(tg <- dbGetQuery(db,statement="select * from gene where Sample in (select Sample from testSamples where Kfold=1) and Gene='MTL5'"))
ts <- dbGetQuery(db,statement="select * from Snps where Sample in (select Sample from testSamples where Kfold=1) and Snps='rs7129419'")

tg

