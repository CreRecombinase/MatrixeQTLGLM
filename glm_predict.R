#Script for creating prediction matrix for eQTLs from glmnet
library(plyr)
library(glmnet)
library(BatchExperiments)
library(reshape2)

dbfile <- "/scratch/nwk2/mEQTL_ERpnc/glmEQTL/unimputed_brca/gene_snp.db"
dbo <- "/scratch/nwk2/mEQTL_ERpnc/glmEQTL/unimputed_brca/predict.db"
db <- dbConnect(drv=dbDriver("SQLite"),dbname=dbfile,loadable.extensions=T)
all.iters <- dbGetQuery(db,"select distinct Gene,Kfold from eqtls")

all.iters <- lapply(chunk(1:nrow(all.iters),n.chunks=200),function(x)all.iters[x,])
glm_predict <- function(t.iters,dbout,dbfile){
  m_ply(.data=t.iters,.fun=function(Kfold,Gene,dbfile,dbo){
    db <- dbConnect(drv=dbDriver("SQLite"),dbname=dbfile)
    querytrain <- paste0("select * from Snps where Sample in (select Sample from trainSamples where Kfold = ",Kfold,") and Snps in (select SNP from eqtls where Gene='",Gene,"' and Kfold=",Kfold,")order by Sample")
    snp.train <- dbGetQuery(db,querytrain)
    snp.train <- acast(data=snp.train,formula=Sample~Snps,value.var="Value")
    querytest <- paste0("select * from Snps where Sample in (select Sample from testSamples where Kfold = ",Kfold,") and Snps in (select SNP from eqtls where Gene='",Gene,"' and Kfold=",Kfold,") order by Sample")
    snp.test <- dbGetQuery(db,querytest)
    snp.test <- acast(data=snp.test,formula=Sample~Snps,value.var="Value")
    exp.query <- paste0("select Value from gene where Sample in (select Sample from trainSamples where Kfold=",Kfold,") and Gene='",Gene,"' order by Sample")
    exp.train <- unlist(dbGetQuery(db,exp.query))
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
        npred <- data.frame(Sample=rownames(tpred),Value=tpred[,1])
        odb <- dbConnect(drv=dbDriver("SQLite"),dbname=dbo)
        db.success <- tryCatch(dbWriteTable(odb,name="predict",npred,row.names=F),error=function(e)e)
        dbDisconnect(odb)
        while(inherits(db.success,"error")){
          odb <- dbConnect(drv=dbDriver("SQLite"),dbname=dbo)
          db.success <- tryCatch(dbWriteTable(odb,name="predict",npred,row.names=F),error=function(e)e)
          dbDisconnect(odb)
        }
      }
    }else{
      tpred <- predict(cv1,newx=snp.test,s=cv1$lambda.1se)
      npred <- data.frame(Sample=rownames(tpred),Value=tpred[,1])
      odb <- dbConnect(drv=dbDriver("SQLite"),dbname=dbo)
      db.success <- tryCatch(dbWriteTable(odb,name="predict",npred,row.names=F,append=T),error=function(e)e)
      dbDisconnect(odb)
      while(inherits(db.success,"error")){
        odb <- dbConnect(drv=dbDriver("SQLite"),dbname=dbo)
        db.success <- tryCatch(dbWriteTable(odb,name="predict",npred,row.names=F,append=T),error=function(e)e)
        dbDisconnect(odb)
      }
    }
    dbDisconnect(db)
    
  },dbfile,dbo,.parallel=F,.paropts=list(.packages=c("RSQLite","glmnet","reshape2"),.inorder=F),.inform=F)
 
    
}


glm.reg <- makeRegistry(registry.name,file.dir=m.dir,packages=c("glmnet","plyr","reshape2"))

batchMap(glm.reg,fun=glm_predict,t.iters=all.iters,more.args=list(dbout=dbo,dbfile=dbfile))



submitJobs(glm.reg)

