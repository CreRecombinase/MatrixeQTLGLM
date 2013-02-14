#Script for creating prediction matrix for eQTLs from glmnet
library(plyr)
library(sqldf)
library(glmnet)
library(BatchExperiments)
library(doParallel)


library(MatrixEQTL)

snp.exploc <- "/scratch/nwk2/mEQTL_ERpnc/glmEQTL/unimputed_brca/snp.exp.Rdata"
sample.num <- 513
kfold <- 57

513/57

train.indices <- chunk(rep(1:sample.num,kfold),n.chunks=kfold)
test.indices <- chunk(1:sample.num,chunk.size=sample.num/kfold)
train.indices <- mapply(FUN=function(x,y)x[-y],train.indices,test.indices,SIMPLIFY=F)

eqtl.base <- "/scratch/nwk2/mEQTL_ERpnc/glmEQTL/unimputed_brca/57-fold/unimputed_brca_trans"
eqtl.files <- paste0(eqtl.base,1:57,".txt")
out.dir <- "/scratch/nwk2/mEQTL_ERpnc/glmEQTL/unimputed_brca/57-fold/"

m.dir <- tempfile("glm_res",tmpdir=out.dir)
registry.name <- paste("glm_reg_")


glm_predict <- function(snp.exploc,train.index,test.index,eqtl.file){
  registerDoParallel(cores=5)
  load(snp.exploc)
  snp.exp$snps <- as.matrix(snp.exp$snps)
  snp.exp$gene <- as.matrix(snp.exp$gene)
  eqtls <- read.csv.sql(eqtl.file,sep="\t",header=T,eol="\n")
  snp.genes <- split(eqtls$SNP,eqtls$gene)
  

  test.train.snp.exp <- foreach(s.g=iter(snp.genes),gn=names(snp.genes),snp.exp$snps,snp.exp$gene,train.index,test.index) %do%{
    snp.train <- t(snp.exp$snps[s.g,train.index,drop=F])
    exp.train <- snp.exp$gene[gn,train.index]
    badcols <- which(is.na(exp.train))
    if(length(badcols)>0){
      snp.train <- snp.train[-badcols,,drop=F]
      exp.train <- exp.train[-badcols]
    }
    snp.test <- t(snp.exp$snps[s.g,test.index,drop=F])
    list(snp.train=snp.train,exp.train=exp.train,snp.test=snp.test,gn=gn)
  }
  checkFun <- function(x){
    return(!all(apply(x$snp.train+1,MARGIN=2,function(x)sum(sort(tabulate(x),decreasing=T)[-1]))<=2))
  }
  predmat <- foreach(s.g=iter(test.train.snp.exp,checkFunc=checkFun),.combine=cbind,.multicombine=T,.inorder=F,.packages="glmnet",.verbose=T,.noexport=c("cv1","fit1")) %dopar% {
    
    cv1 <- cv.glmnet(x=s.g$snp.train,s.g$exp.train);
    fit1 <- glmnet(s.g$snp.train,s.g$exp.train,lambda=cv1$lambda.1se,alpha=0.95);
    tpred <- predict(fit1,newx=s.g$snp.test);
    colnames(tpred)<- s.g$gn;
    message(paste0(Sys.time(),"\n"))
    tpred
    
  }
  return(predmat)
}




glm.reg <- makeRegistry(registry.name,file.dir=m.dir,packages=c("doParallel","glmnet","sqldf"))

batchMap(glm.reg,glm_predict,train.index=train.indices,test.index=test.indices,eqtl.file=eqtl.files,more.args=list(snp.exploc=snp.exploc))



submitJobs(glm.reg)

