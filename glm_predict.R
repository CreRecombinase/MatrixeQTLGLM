#Script for creating prediction matrix for eQTLs from glmnet
## library(plyr,quietly=T)
library(glmnet,quietly=T)
library(BatchExperiments,quietly=T)
library(reshape2,quietly=T)
library(RSQLite,quietly=T)
library(doParallel,quietly=T)

#usage glm_predict.R <DBFILE> <chunks> <kfold> <out.dir> <queue> <memory> <time> <threads>
makeClusterFunctionsLSF("~/lsf-threaded.tmpl")

oargs <- commandArgs(trailingOnly=TRUE)


dbfile <- oargs[1]
chunks <- as.integer(oargs[2])
max.kfolds <- as.integer(oargs[3])
out.dir <- oargs[4]
queue <- oargs[5]
memory <- as.integer(oargs[6])
time <- oargs[7]
threads <- as.integer(oargs[8])



db <- dbConnect(drv=dbDriver("SQLite"),dbname=dbfile,loadable.extensions=T)


all.iters <- dbGetQuery(db,"select distinct gene from ggenes")[[1]]

dbDisconnect(db)
all.iters <- chunk(all.iters,n.chunks=chunks)


glm_predict <- function(ot.iters,dbfile,threads,kfolds){
  
  registerDoParallel(cores=threads-1)
  
  create.set <- function(items){
    items <-paste0("'",items,"'")
    paste0("(",Reduce(function(x,y)paste(x,y,sep=","),items),")")
  }
  
  db <- dbConnect(drv=dbDriver("SQLite"),dbname=dbfile,loadable.extensions=T)
  snpsql <- "select SNP,gene,Kfold from eqtls where gene = :gene"
  snps <- dbGetPreparedQuery(db,snpsql,data.frame(gene=ot.iters))
  all.exp <- acast(dbGetQuery(db,paste0("select * from gene where Gene in ",create.set(ot.iters))),formula=Sample~Gene,value.var="Value")
  nnsql <- paste0("select * from snps where Snp in ",create.set(unique(snps$SNP)))
  nnsnps <- dbGetQuery(db,nnsql)
  trainSamples <- dbGetQuery(db,"select Sample,Kfold from trainSamples")
  trainSamples <- split(trainSamples$Sample,trainSamples$Kfold)
  testSamples <- dbGetQuery(db,"select Sample,Kfold from testSamples")
  testSamples <- split(testSamples$Sample,testSamples$Kfold)
  asnps <- acast(data=nnsnps,formula=Sample~Snp,value.var="Value")
  eqtls <- sapply(split(snps,snps$gene),FUN=function(x)split(x$SNP,x$Kfold),simplify=F)
  
  
  
  glm.engine <- function(gene,Kfold){
    snp.train <- asnps[trainSamples[[Kfold]],eqtls[[gene]][[Kfold]],drop=F]
    snp.test<- asnps[testSamples[[Kfold]],eqtls[[gene]][[Kfold]],drop=F]
    exp.train <- all.exp[trainSamples[[Kfold]],gene]
    cv1 <- tryCatch(cv.glmnet(x=snp.train,exp.train,alpha=0.95),error=function(e)e)
    if(inherits(cv1,"error")){
      if(any(apply(snp.train,MARGIN=2,function(x)sum(sort(tabulate(x),decreasing=T)[-1]))>=2)){
        badsamples <- is.na(exp.train)
        snp.train <- snp.train[!badsamples,,drop=F]
        exp.train <- exp.train[!badsamples]
        t <- 0
        while(inherits(cv1,"error")&&t<3){
          cv1 <- tryCatch(cv.glmnet(x=snp.train,exp.train,alpha=0.95),error=function(e)e)
          t <- t+1
        }
        if(t>=3){
          return(NULL)
        }
        tpred <- predict(cv1,newx=snp.test,s=cv1$lambda.1se)
        npred <- data.frame(Sample=rownames(tpred),Value=tpred[,1],Gene=gene)
        cf <- as.matrix(coef(cv1,s=cv1$lambda.1se))
        return(list(pred=npred,coefs=cf,lambda=cv1$lambda.1se))
      }     
    }else{
      tpred <- predict(cv1,newx=snp.test,s=cv1$lambda.1se)
      npred <- data.frame(Sample=rownames(tpred),Value=tpred[,1],Gene=gene)
      cf <- as.matrix(coef(cv1,s=cv1$lambda.1se))
      return(list(pred=npred,coefs=cf,lambda=cv1$lambda.1se))
     
    }
    
  }

  t.iters <- expand.grid(gene=ot.iters,Kfold=as.character(1:kfolds),stringsAsFactors=F)
  
  
  tt.res <- mlply(.data=t.iters,.fun=glm.engine,.parallel=T,.inform=F,.paropts=list(.multicombine=T,.inorder=F,.verbose=F,.export=c("glm.engine","asnps","all.exp","eqtls","trainSamples","testSamples"),.packages=c("glmnet")))
  return(tt.res)
    
}





m.dir <- tempfile("glm_res",tmpdir=out.dir)

glm.reg <- makeRegistry("glmreg",file.dir=m.dir,packages=c("glmnet","plyr","reshape2","RSQLite","doParallel"))

batchMap(glm.reg,fun=glm_predict,ot.iters=all.iters,more.args=list(dbfile=dbfile,threads=threads,kfolds=max.kfolds))


## submitJobs(glm.reg,resources=list(queue=queue,threads=threads,memory=memory,time=time))

## Sys.sleep(10)
