#Script for creating prediction matrix for eQTLs from glmnet
library(plyr)
library(glmnet)
library(BatchExperiments)
library(reshape2)
library(RSQLite)


if(Sys.info()['sysname']=="Windows"){
  root.dir <- "C:/Users/nknoblau/Documents/R_WS/MatrixeQTLGLM/"
  out.dir <- root.dir
  dbfile <- "D:/gene_snp.db"
  dbo <- "D:/predict_gene.db"
  chunks <- 10
}else{
  root.dir <- "/scratch/nwk2/mEQTL_ERpnc/glmEQTL/unimputed_brca/"
  out.dir <- paste0(root.dir,"57-fold")
  dbfile <- "/scratch/nwk2/mEQTL_ERpnc/glmEQTL/unimputed_brca/gene_snp.db"
  dbo <- "/scratch/nwk2/mEQTL_ERpnc/glmEQTL/unimputed_brca/predict.db"
  chunks=200
}



db <- dbConnect(drv=dbDriver("SQLite"),dbname=dbfile,loadable.extensions=T)
all.iters <- dbGetQuery(db,"select distinct Gene,Kfold from eqtls")

if(Sys.info()['sysname']=="Windows"){
  all.iters <- all.iters[1:100,]
}
  
all.iters <- lapply(chunk(1:nrow(all.iters),n.chunks=chunks),function(x)all.iters[x,])


glm_predict <- function(t.iters,dbo,dbfile){
  
  glm.engine <- function(Kfold,Gene){
    db <- dbConnect(drv=dbDriver("SQLite"),dbname=dbfile)
    querytrain <- paste0("select * from Snps where Sample in (select Sample from trainSamples where Kfold = ",Kfold,") and Snps in (select SNP from eqtls where Gene='",Gene,"' and Kfold=",Kfold,")order by Sample")
    snp.train <- dbGetQuery(db,querytrain)
    snp.train <- acast(data=snp.train,formula=Sample~Snps,value.var="Value")
    querytest <- paste0("select * from Snps where Sample in (select Sample from testSamples where Kfold = ",Kfold,") and Snps in (select SNP from eqtls where Gene='",Gene,"' and Kfold=",Kfold,") order by Sample")
    snp.test <- dbGetQuery(db,querytest)
    snp.test <- acast(data=snp.test,formula=Sample~Snps,value.var="Value")
    exp.query <- paste0("select Value from gene where Sample in (select Sample from trainSamples where Kfold=",Kfold,") and Gene='",Gene,"' order by Sample")
    exp.train <- unlist(dbGetQuery(db,exp.query))
    dbDisconnect(db)
    cv1 <- tryCatch(cv.glmnet(x=snp.train,exp.train,alpha=0.95),error=function(e)e)
    if(inherits(cv1,"error")){
      if(sum(sort(tabulate(snp.train),decreasing=T)[-1])>=2){
        t <- 0
        while(inherits(cv1,"error")||t<10){
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
  
  
  tt.res <- mdply(.data=t.iters,.fun=glm.engine)
  tt.res <- tt.res[,c("Sample","Value","Gene")]
  return(tt.res)
    
}





m.dir <- tempfile("glm.res",tmpdir=out.dir)

glm.reg <- makeRegistry("glmreg",file.dir=m.dir,packages=c("glmnet","plyr","reshape2","RSQLite"))

batchMap(glm.reg,fun=glm_predict,t.iters=all.iters,more.args=list(dbo=dbo,dbfile=dbfile))


submitJobs(glm.reg)

