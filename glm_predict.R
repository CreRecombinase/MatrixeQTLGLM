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



train.indices <- chunk(rep(1:sample.num,kfold),n.chunks=kfold)
test.indices <- chunk(1:sample.num,chunk.size=sample.num/kfold)
train.indices <- mapply(FUN=function(x,y)x[-y],train.indices,test.indices,SIMPLIFY=F)

eqtl.base <- "/scratch/nwk2/mEQTL_ERpnc/glmEQTL/unimputed_brca/57-fold/unimputed_brca_trans"
eqtl.files <- paste0(eqtl.base,1:57,".txt")
out.dir <- "/scratch/nwk2/mEQTL_ERpnc/glmEQTL/unimputed_brca/57-fold/"

m.dir <- tempfile("glm_res",tmpdir=out.dir)
registry.name <- paste("glm_reg_")


glm_predict <- function(snp.exploc,train.index,test.index,eqtl.file){
  registerDoParallel(cores=12)
  load(snp.exploc)
  snp.exp$snps <- as.matrix(snp.exp$snps)
  snp.exp$gene <- as.matrix(snp.exp$gene)
  eqtls <- read.csv.sql(eqtl.file,sep="\t",header=T,eol="\n")
  snp.genes <- split(eqtls$SNP,eqtls$gene)
  
  checkFun <- function(x){
    return(!all(apply(x,MARGIN=2,function(x)sum(sort(tabulate(x),decreasing=T)[-1]))<=2))
  }
  
  test.train.snp.exp <- lapply(1:length(snp.genes),function(i){
    snp.train <- t(snp.exp$snps[snp.genes[[i]],train.index,drop=F])
    exp.train <- snp.exp$gene[names(snp.genes)[i],match(rownames(snp.train),colnames(snp.exp$gene))]
    badcols <- which(is.na(exp.train))
    if(length(badcols)>0){
      snp.train <- snp.train[-badcols,,drop=F]
      exp.train <- exp.train[-badcols]
    }
    snp.test <- t(snp.exp$snps[snp.genes[[i]],test.index,drop=F])
    if(checkFun(snp.train+1)){
      return(list(snp.train=snp.train,exp.train=exp.train,snp.test=snp.test,gn=names(snp.genes)[i]))
    }
  })
  test.train.snp.exp <- test.train.snp.exp[!sapply(test.train.snp.exp,is.null)]

 
  genpred <- function(s.g){
    cv1 <- cv.glmnet(x=s.g$snp.train,s.g$exp.train,alpha=0.95)
    tpred <- predict(cv1,newx=s.g$snp.test,s=cv1$lambda.1se)
    colnames(tpred)<-s.g$gn
    return(data.frame(t(tpred)))
  }
  
  predmat <- ldply(test.train.snp.exp,.fun=genpred,.parallel=T,.paropts=list(.packages="glmnet"))
  return(predmat)
}


glm.reg <- makeRegistry(registry.name,file.dir=m.dir,packages=c("doParallel","glmnet","sqldf","plyr"))

batchMap(glm.reg,glm_predict,train.index=train.indices,test.index=test.indices,eqtl.file=eqtl.files,more.args=list(snp.exploc=snp.exploc))



submitJobs(glm.reg)

