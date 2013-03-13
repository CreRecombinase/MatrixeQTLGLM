#Code to generate list of lists to be fed into glm_predict.R

library(MatrixEQTL)
library(plyr)
library(doParallel)
library(sqldf)
library(RSQLite)
library(RSQLite.extfuns)
library(BatchExperiments)
library(reshape2)



base.dir <- "/scratch/nwk2/mEQTL_ERpnc/glmEQTL/brca_RNAseq/"

dbfile <- paste0(base.dir,"rnaseq_snp_genes.db")

eqtl.files <- paste0(base.dir,"82-fold/unimputed_brca_RNAseq_trans")
eqtl.files <- paste0(eqtl.files,1:82,".txt")

exp.file <- paste0(base.dir,"brca_RNAseq_expression.txt")
snp.file <- paste0(base.dir,"unimputed_brca_RNAseq_snp.txt")
sample.num <- 819
kfold <- 82



train.indices <- chunk(rep(1:sample.num,kfold),n.chunks=kfold)
test.indices <- chunk(1:sample.num,chunk.size=ceiling(sample.num/kfold))
train.indices <- mapply(FUN=function(x,y)x[-y],train.indices,test.indices,SIMPLIFY=F)




aeqtls <- ldply(eqtl.files,function(x){
  tqtl <- read.csv.sql(x,header=T,eol="\n",sep="\t")
  knum <- gsub(".+trans([0-9]+).txt","\\1",x)
  tqtl$knum <- knum
  return(tqtl)
}
)

db <- dbConnect(drv=dbDriver("SQLite"),dbname=dbfile,loadable.extensions=T)
colnames(aeqtls) <- c("SNP","Gene","tStat","pValue","FDR","Kfold")
dbWriteTable(db,"eqtls",aeqtls,row.names=F)
dbSendQuery(db,"Create index sgk on eqtls(SNP,Gene,Kfold)")
dbSendQuery(db,"Create index gk on eqtls(Gene,Kfold)")



fsnp.genes <- split(aeqtls,aeqtls$knum)
snp.genes <- unlist(lapply(fsnp.genes,function(x)split(x$SNP,x$gene)),recursive=F)
rm(aeqtls,fsnp.genes)

snp.cases <- scan(snp.file,what="character",nlines=1,sep="\n")
snp.cases <- strsplit(snp.cases,split="\t")[[1]]
snp.cases <- snp.cases[-1]
exp.cases <- scan(exp.file,what="character",nlines=1,sep="\n")
exp.cases <- strsplit(exp.cases,"\t")[[1]]
exp.cases <- exp.cases[-1]
if(any(exp.cases!=snp.cases)){
  stop("Problem with MatrixEQTL inputs")
  
}


test.df <- do.call("rbind",lapply(1:length(test.indices),function(x){
  data.frame(Sample=exp.cases[test.indices[[x]]],Kfold=x,Index=test.indices[[x]])
}))

train.df <-test.df <- do.call("rbind",lapply(1:length(train.indices),function(x){
  data.frame(Sample=exp.cases[train.indices[[x]]],Kfold=x,Index=train.indices[[x]])
}))

dbWriteTable(db,name="testSamples",test.df,row.names=F,overwrite=T,append=F)
dbWriteTable(db,name="trainSamples",train.df,row.names=F,overwrite=T,append=F)
dbSendQuery(db,"Create index k on testSamples(Kfold)")
dbSendQuery(db,"Create index ktr on trainSamples(Kfold)")

dbDisconnect(db)



