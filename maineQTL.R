#First revision of Andy's pseudocode for unimputed BRCA SNP and Expression Data
#2/6/13
#NWK
library(sqldf)
library(plyr)
library(BatchExperiments)
library(MatrixEQTL)
x <- SlicedData$

root.dir <- "C:/Users/nknoblau/Documents/R_WS/MatrixEQTL_ERpcn/"
setwd(root.dir)
snp.filepath <- "SNP6.txt"
exp.filepath <- "nLevel3.txt"

snp.loc.fp <- "matrixSNPlocations.txt"
exp.loc.fp <- "L3genelocs.txt"


brca.dat <- list(snp=SlicedData$new(load.data.matrix(snp.filepath)),exp=SlicedData$new(load.data.matrix(exp.filepath)),exp.loc=load.anno(exp.loc.fp),snp.loc=load.anno(snp.loc.fp))







s1<-sample(c(1:n),pts,replace=T)
for(i in 1:n){
  exp.train<-exp[s1!=i,]
  snp.train<-snp[s1!=i,]
  exp.test<-exp[s1==i,]
  snp.test<-snp[s1==i,]
  
  matrix.eqtl.out<-matrix.eqtl(snp.train,exp.train)
  activeSNPs<-SNPs that are in at least one interaction in matrix.eqtl.out at some threshold
  activeGenes<-genes in at least one interaction
  for(j in 1:length(activeGenes)){
    x=snp.train[,activeSNPsForGenej]
    y=exp.train[,activeGenes[j]]
    cv1<-cv.glmnet(x=x,y=y)
    fit1<-glmnet(x=x,y=y,lambda=cv1$lambda.1se)
    test.x<-snp.test[,activeSNPsForGenej]
    predMatrix[s1==i,colnames(predMatrix)[j]]<-predict(fit1,newx=test.x)
  }
}
## Now, we can compute mean squared errors...   






load.data.matrix <- function(filepath){
  rawdat <- read.csv.sql(filepath,sep="\t",header=T,eol="\n")
  rownames(rawdat)<- rawdat[,1]
  rawdat <- rawdat[,-1]
  rawdat <- data.matrix(rawdat)
  colnames(rawdat)<- gsub("_","-",colnames(rawdat))
  return(rawdat)
}
load.anno <- function(filepath){
  read.csv.sql(filepath,sep="\t",header=T,eol="\n")
}


