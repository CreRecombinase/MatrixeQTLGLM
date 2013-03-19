#Script for creating prediction matrix for eQTLs from glmnet
library(plyr,quietly=T)
library(glmnet,quietly=T)
library(BatchExperiments,quietly=T)
library(reshape2,quietly=T)
library(RSQLite,quietly=T)
library(compiler,quietly=T)
library(doParallel,quietly=T)


if(Sys.info()['sysname']=="Windows"){
  root.dir <- "C:/Users/nknoblau/Documents/R_WS/MatrixeQTLGLM/testdata/"
  out.dir <- root.dir
  dbfile <- "C:/Users/nknoblau/Documents/R_WS/MatrixeQTLGLM/testdata/testdb.db"
  chunks <- 10
}else{
  root.dir <- "/scratch/nwk2/mEQTL_ERpnc/glmEQTL/unimputed_brca/"
  out.dir <- paste(root.dir,"57-fold",sep="")
  dbfile <- "/scratch/nwk2/mEQTL_ERpnc/glmEQTL/unimputed_brca/gene_snp.db"
  chunks=100
}



db <- dbConnect(drv=dbDriver("SQLite"),dbname=dbfile,loadable.extensions=T)

all.iters <- dbGetQuery(db,"select distinct eqtls.gene,Kfold from eqtls,ggenes where eqtls.gene=ggenes.Gene")

dbDisconnect(db)


  
all.iters <- lapply(chunk(1:nrow(all.iters),n.chunks=chunks),function(x)all.iters[x,])


glm_predict <- function(t.iters,dbfile){

  registerDoParallel(cores=4)
  
  glm.engine <- function(gene,Kfold){
    
    db <- dbConnect(drv=dbDriver("SQLite"),dbname=dbfile,loadable.extensions=T)
    dbGetQuery(db,"pragma journal_mode=TRUNCATE")
    dbGetQuery(db,"pragma busy_timeout=20000")
    querytrain <- paste("select * from snps where Sample in (select Sample from trainSamples where Kfold = ",Kfold,") and Snp in (select SNP from eqtls where Gene='",gene,"' and Kfold=",Kfold,")order by Sample",sep="")
    snp.train <- dbGetQuery(db,querytrain)
    snp.train <- acast(data=snp.train,formula=Sample~Snp,value.var="Value")
    querytest <- paste("select * from snps where Sample in (select Sample from testSamples where Kfold = ",Kfold,") and Snp in (select SNP from eqtls where Gene='",gene,"' and Kfold=",Kfold,") order by Sample",sep="")
    snp.test <- dbGetQuery(db,querytest)
    snp.test <- acast(data=snp.test,formula=Sample~Snp,value.var="Value")
    exp.query <- paste("select Value from gene where Sample in (select Sample from trainSamples where Kfold=",Kfold,") and Gene='",gene,"' order by Sample",sep="")
    exp.train <- unlist(dbGetQuery(db,exp.query))
    dbDisconnect(db)
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
        if(t>=10){
          return(NULL)
        }
        tpred <- predict(cv1,newx=snp.test,s=cv1$lambda.1se)
        npred <- data.frame(Sample=rownames(tpred),Value=tpred[,1],Gene=gene)
        return(npred)
      }     
    }else{
      tpred <- predict(cv1,newx=snp.test,s=cv1$lambda.1se)
      npred <- data.frame(Sample=rownames(tpred),Value=tpred[,1],Gene=gene)
      return(npred)
     
    }
    
  }
  
  nglm.engine <- cmpfun(glm.engine,options=list(optimize=3))
  system.time(tt.res <- mdply(.data=t.iters,.fun=nglm.engine,.parallel=T,.paropts=list(.multicombine=T,.inorder=F,.verbose=F,.export=c("nglm.engine","dbfile"),.packages=c("glmnet","RSQLite","reshape2"))))
  tt.res <- tt.res[,c("Sample","Value","Gene")]
  return(tt.res)
    
}





m.dir <- tempfile("glm_res",tmpdir=out.dir)

glm.reg <- makeRegistry("glmreg",file.dir=m.dir,packages=c("glmnet","plyr","reshape2","RSQLite","doParallel"))

batchMap(glm.reg,fun=glm_predict,t.iters=all.iters,more.args=list(dbfile=dbfile))


submitJobs(glm.reg)


