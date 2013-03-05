#Script for creating prediction matrix for eQTLs from glmnet
library(plyr)
library(sqldf)
library(glmnet)
library(BatchExperiments)
library(doParallel)
library(reshape2)

dbfile <- "D:/gene_snp.db"
db <- dbConnect(drv=dbDriver("SQLite"),dbname=dbfile,loadable.extensions=T)
agenes <- dbGetQuery(db,"select distinct Gene,Kfold from eqtls")





all.iters <- lapply(chunk(1:nrow(agenes),n.chunks=2000),function(x)agenes[x,])

t.iters <- all.iters[[1]]
colnames(t.iters)<- c("gn","kf")
gn <- t.iters[1,1]
kf <- as.integer(t.iters[1,2])
glm_predict <- function(db,train.index,test.index){
  registerDoParallel(cores=6)

  system.time(a.resultes <- mlply(t.iters[1:50,],.fun=function(kf,gn,dbfile){
    db <- dbConnect(drv=dbDriver("SQLite"),dbname=dbfile)
    querytrain <- paste0("select * from Snps where Sample in (select Sample from trainSamples where Kfold = ",kf,") and Snps in (select SNP from eqtls where Gene='",gn,"' and Kfold=",kf,")order by Sample")
    snp.train <- dbGetQuery(db,querytrain)
    snp.train <- acast(data=snp.train,formula=Sample~Snps)
    querytest <- paste0("select * from Snps where Sample in (select Sample from testSamples where Kfold = ",kf,") and Snps in (select SNP from eqtls where Gene='",gn,"' and Kfold=",kf,") order by Sample")
    snp.test <- dbGetQuery(db,querytest)
    snp.test <- acast(data=snp.test,formula=Sample~Snps)
    exp.query <- paste0("select Value from gene where Sample in (select Sample from trainSamples where Kfold=",kf,") and Gene='",gn,"' order by Sample")
    exp.train <- unlist(dbGetQuery(db,exp.query))
    cv1 <- cv.glmnet(x=snp.train,exp.train,alpha=0.95)
    tpred <- data.frame(predict(cv1,newx=snp.test,s=cv1$lambda.1se))
    colnames(tpred)<-gn
    dbDisconnect(db)
    return(tpred)
  },dbfile,.parallel=T,.paropts=list(.packages=c("RSQLite","glmnet","reshape2")),.inform=F))
  
  nbsnps <- dbGetQuery(db,"select * from Snps where Snps in (select SNP from eqtls where Gene='15E1.2' and Kfold=39)")
  test.cases <-  dbGetQuery(db,"select distinct Sample from testSamples")
  train.cases <- dbGetQuery(db,"select distinct Sample from trainSamples")
  
  
  
  
  
  
 
  genpred <- function(s.g){
    cv1 <- cv.glmnet(x=s.g$snp.train,s.g$exp.train,alpha=0.95)
    tpred <- predict(cv1,newx=s.g$snp.test,s=cv1$lambda.1se)
    tpred <- data.frame(t(tpred))
    tpred$gn <- s.g$gn
    return(tpred)
  }
  
  predmat <- ldply(test.train.snp.exp,.fun=genpred,.parallel=T,.paropts=list(.packages="glmnet",.multicombine=T,.inorder=F))
  
  system.time(tpredmat <- ldply(test.train.snp.exp[1:300],.fun=genpred,.parallel=T,.paropts=list(.packages="glmnet",.multicombine=T,.inorder=F)))
  Rprof(NULL)
  summaryRprof("glmprof.out")
  tpredmat
  return(predmat)
}


glm.reg <- makeRegistry(registry.name,file.dir=m.dir,packages=c("doParallel","glmnet","sqldf","plyr","MatrixEQTL"))

batchMap(glm.reg,glm_predict,train.index=train.indices,test.index=test.indices,eqtl.file=eqtl.files,more.args=list(snp.exploc=snp.exploc))



submitJobs(glm.reg)

