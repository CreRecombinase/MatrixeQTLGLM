#Script for creating prediction matrix for eQTLs from glmnet
library(plyr)
library(sqldf)
library(glmnet)
library(BatchExperiments)
library(doParallel)




glm_predict <- function(snp.exploc,train.index,test.index){
  registerDoParallel(cores=12)


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

