#Script for creating prediction matrix for eQTLs from glmnet
library(plyr)
library(glmnet)
library(BatchExperiments)
library(reshape2)
library(RSQLite)
library(doParallel)


if(Sys.info()['sysname']=="Windows"){
  root.dir <- "C:/Users/nknoblau/Documents/R_WS/MatrixeQTLGLM/"
  out.dir <- root.dir
  dbfile <- "D:/gene_snp.db"
  chunks <- 100
}else{
  root.dir <- "/scratch/nwk2/mEQTL_ERpnc/glmEQTL/unimputed_brca/"
  out.dir <- paste(root.dir,"57-fold",sep="")
  dbfile <- "/scratch/nwk2/mEQTL_ERpnc/glmEQTL/unimputed_brca/gene_snp.db"
  chunks=100
}



db <- dbConnect(drv=dbDriver("SQLite"),dbname=dbfile,loadable.extensions=T)

all.iters <- dbGetQuery(db,"select distinct Gene,Kfold from eqtls")
dbDisconnect(db)


  
all.iters <- lapply(chunk(1:nrow(all.iters),n.chunks=chunks),function(x)all.iters[x,])


glm_predict <- function(t.iters,dbfile){

  registerDoParallel(cores=4)
  
  glm.engine <- function(Kfold,Gene){
    
    db <- dbConnect(drv=dbDriver("SQLite"),dbname=dbfile,loadable.extensions=T)
    dbGetQuery(db,"pragma journal_mode=TRUNCATE")
    querytrain <- paste("select * from Snps where Sample in (select Sample from trainSamples where Kfold = ",Kfold,") and Snps in (select SNP from eqtls where Gene='",Gene,"' and Kfold=",Kfold,")order by Sample",sep="")
    snp.train <- dbGetQuery(db,querytrain)
    snp.train <- acast(data=snp.train,formula=Sample~Snps,value.var="Value")
    querytest <- paste("select * from Snps where Sample in (select Sample from testSamples where Kfold = ",Kfold,") and Snps in (select SNP from eqtls where Gene='",Gene,"' and Kfold=",Kfold,") order by Sample",sep="")
    snp.test <- dbGetQuery(db,querytest)
    snp.test <- acast(data=snp.test,formula=Sample~Snps,value.var="Value")
    exp.query <- paste("select Value from gene where Sample in (select Sample from trainSamples where Kfold=",Kfold,") and Gene='",Gene,"' order by Sample",sep="")
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
        npred <- data.frame(Sample=rownames(tpred),Value=tpred[,1],Gene=Gene)
        return(npred)
      }     
    }else{
      tpred <- predict(cv1,newx=snp.test,s=cv1$lambda.1se)
      npred <- data.frame(Sample=rownames(tpred),Value=tpred[,1],Gene=Gene)
      return(npred)
     
    }
    
  }
  
  
  system.time(tt.res <- mdply(.data=t.iters,.fun=glm.engine,.parallel=T,.paropts=list(.multicombine=T,.inorder=F,.verbose=F,.export=c("glm.engine","dbfile"),.packages=c("glmnet","RSQLite","reshape2"))))
  tt.res <- tt.res[,c("Sample","Value","Gene")]
  return(tt.res)
    
}





#m.dir <- tempfile("glm_res",tmpdir=out.dir)

#glm.reg <- makeRegistry("glmreg",file.dir=m.dir,packages=c("glmnet","plyr","reshape2","RSQLite","doParallel"))

#batchMap(glm.reg,fun=glm_predict,t.iters=all.iters,more.args=list(dbfile=dbfile))


#submitJobs(glm.reg)


