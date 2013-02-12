#Script for creating prediction matrix for eQTLs from glmnet
library(plyr)
library(sqldf)
library(glmnet)
library(BatchExperiments)

library(MatrixEQTL)
snp.exploc <- "C:/Users/nknoblau/Documents/R_WS/MatrixeQTLGLM/test/snp.exp.Rdata"
load(snp.exploc)
snp.exp$snps <- as.matrix(snp.exp$snps)

glm_predict <- function(outfiles,snp.exploc,train.index,test.index,eqtl.file){
  load(snp.expfile)
  snp.exp$snps <- as.matrix(snp.exp$snps)
  snp.exp$gene <- as.matrix(snp.exp$gene)
  eqtls <- read.csv.sql(eqtl.file,sep="\t",header=T,eol="\n")
  snp.genes <- split(eqtls$SNP,eqtls$gene)
  
  
  
  
  
  
  
  
  
  

  ngenes <- nrow(snp.exp$gene)
  nSamples <- ncol(snp.exp$gene)
  pred.Mat <- matrix(NA,nrow=nSamples,ncol=ngenes,dimnames=list(colnames(snp.exp$snps),rownames(snp.exp$gene)))
  eqtl.dir <- "C:/Users/nknoblau/Documents/R_WS/MatrixeQTLGLM/test/"
  eqtl.res <- dir("C:/Users/nknoblau/Documents/R_WS/MatrixeQTLGLM/test/",pattern=".+trans.+")
  all.tra.eqtls <- ldply(eqtl.res,.fun=function(x){
    tres <- read.csv.sql(paste(eqtl.dir,x,sep=""),sep="\t",header=T,eol="\n")
    tres$fn <- x
    return(tres)
  })
  all.tra.eqtls$fn <- gsub(".+([0-9]).+","\\1",all.tra.eqtls$fn)
  all.tra.eqtls$fn <- as.integer(all.tra.eqtls$fn)
  
}
tsnp1 <- split(x=all.tra.eqtls[all.tra.eqtls$fn=="1","SNP"],all.tra.eqtls[all.tra.eqtls$fn=="1","gene"])





make.pred <- function(train.snps,train.gene,test.snps,train.ind,test.ind,snpdat,exp){
  if(ncol(train.snps)==1){
    if(sum(sort(table(train.snps[,1]),decreasing=T)[-1])<=2){
      return(matrix(NA,nrow=nrow(test.snps),ncol=1,dimnames=list(rownames(test.snps),"s0")))
    }
  }
  tables <- apply(train.snps,2,function(x)sort(table(x),decreasing=T))
  table.sums <- sapply(tables,function(x)sum(x[-1]))
  if(all(table.sums<=2)){
    return(matrix(NA,nrow=nrow(test.snps),ncol=1,dimnames=list(rownames(test.snps),"s0")))
  }
  badcols <- which(is.na(train.gene))
  if(length(badcols)>0){
    train.snps <- train.snps[-badcols,,drop=F]
    train.gene <- train.gene[-badcols]
  }
  cv1 <- cv.glmnet(x=train.snps,y=train.gene)
  fit1 <- glmnet(train.snps,train.gene,lambda=cv1$lambda.1se,alpha=0.95)
  tpred <- predict(fit1,newx=test.snps)
  
  return(tpred)
}

zero.var <- function(mat){
  if(ncol(mat)==1){
    return(sum(sort(table(mat[,1]),decreasing=T)[-1])<=2)
  }
  tables <- apply(mat,2,function(x)sort(table(x),decreasing=T))
  table.sums <- sapply(tables,function(x)sum(x[-1]))
  return(all(table.sums<=2))
}


i <- 1
kqtls <- split(all.tra.eqtls,all.tra.eqtls$fn)
dlsnpl <- sapply(kqtls,function(x)split(x$SNP,x$gene))

for(i in 3:length(unique(all.tra.eqtls$fn))){
  tsnpl <- split(x=all.tra.eqtls[all.tra.eqtls$fn==i,"SNP"],all.tra.eqtls[all.tra.eqtls$fn==i,"gene"])
  train.ind <- train.indices[[i]]
  test.ind  <- test.indices[[i]]
  for( j in 1:length(tsnpl)){
    res.snps <- tsnpl[[j]]
    res.gene <- names(tsnpl)[j]
    train.snps <- t(data.matrix(snp.exp$snps[res.snps,train.ind,drop=F]))
    train.gene <- as.vector(snp.exp$gene[res.gene,train.ind])
    badcols <- which(is.na(train.gene))
    if(length(badcols)>0){
      train.snps <- train.snps[-badcols,,drop=F]
      train.gene <- train.gene[-badcols]
    }
    #Look out for "zero variance predictors"
    if(!zero.var(train.snps)){
      test.snps <- t(data.matrix(snp.exp$snps[res.snps,test.ind,drop=F]))
      tpred <- tryCatch(make.pred(train.snps,train.gene,test.snps,train.ind,test.ind),error=function(e)e)
      t <- 0
      while(inherits(tpred,"error")){
        tpred <- tryCatch(make.pred(train.snps,train.gene,test.snps,train.ind,test.ind),error=function(e)e)
        t <- t+1
        if(t==10){
          break
        }
      }
      if(t!=10){
        pred.Mat[rownames(tpred),res.gene]<-tpred
      }
      
    }
    if(j%%100==0){
      print(paste0("j=",j))
    }
  }
  print(paste0("I=",i))
}


bad.rows <- apply(pred.Mat,2,function(x)all(is.na(x)))
sum(bad.rows)
bad.cols <- apply(pred.Mat,1,function(x)all(is.na(x)))
sum(bad.cols)
good.pred <- pred.Mat[,!bad.rows]

head(good.pred)
