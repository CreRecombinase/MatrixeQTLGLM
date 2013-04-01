library(RSQLite)
library(BatchJobs)
library(plyr)

#usage import_results.R <results.dir> <database.file>


oargs <- commandArgs(trailingOnly=T)
results.dir <- oargs[1]
database.file <- oargs[2]
db <- dbConnect(drv=dbDriver("SQLite"),dbname=database.file)

glm.reg <- loadRegistry(results.dir)
testSamples <- dbReadTable(db,"testSamples")
rownames(testSamples) <- testSamples$Sample

for(i in 1:nrow(getJobInfo(glm.reg))){
  tx <- loadResult(glm.reg,i)
  tx <- tx[which(sapply(tx,function(x)!is.null(x)))]
  xdf <- do.call("rbind",sapply(tx,function(x)x$pred,simplify=F))
  dbWriteTable(db,"prediction",xdf,row.names=F,append=T,overwrite=F)
  genes <- sapply(tx,function(x)as.character(unique(x$pred$Gene)))
  kfolds <- testSamples[sapply(tx,function(x)x$pred[1,1]),"Kfold"]
  acoefs <- sapply(tx,function(x)x$coefs[x$coefs!=0,,drop=F])
  bcoef <- do.call("rbind",
                   mapply(FUN=function(coefs,kfolds,genes)
                          {
                            return(data.frame(SNP=rownames(coefs),Gene=genes,Kfolds=kfolds,Value=coefs[,1],stringsAsFactors=F))
                          },
                          coefs=acoefs,genes=genes,kfolds=kfolds,SIMPLIFY=F)
                   )
  lambdas <- data.frame(Gene=genes,Kfold=kfolds,Value=sapply(tx,function(x)x$lambda),stringsAsFactors=F)
  cvm <- data.frame(Gene=genes,Kfold=kfolds,mincvm=sapply(tx,function(x)x$min.cvm),maxcvm=sapply(tx,function(x)x$max.cvm))
  dbWriteTable(db,"cvm",cvm,append=T,overwrite=F,row.names=F)
  dbWriteTable(db,"coefficients",bcoef,append=T,overwrite=F,row.names=F)
  dbWriteTable(db,"lambda",lambdas,append=T,row.names=F,overwrite=F)
  print(i)
}
  
  
