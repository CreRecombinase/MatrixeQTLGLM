#First revision of Andy's pseudocode for unimputed BRCA SNP and Expression Data
#2/6/13
#NWK
library(sqldf)
library(plyr)
library(BatchExperiments)
library(MatrixEQTL)
x <- SlicedData$

root.dir <- "/scratch/nwk2/mEQTL_ERpnc/glmEQTL/unimputed_brca"
setwd(root.dir)
snp.filepath <- "unimputed_brca_snp.txt"
exp.filepath <- "brca_expression.txt"

snp.loc.fp <- "unimputed_brca_snp_anno.txt"
exp.loc.fp <- "brca_expression_anno.txt"

datfile <- "static.Rdata"

if(!file.exists(datfile)){
  shell(paste("/home/nwk2/glm_eqtl/MatrixeQTLGLM/load_static.R","F",snp.filepath,exp.filepath,snp.loc.fp,exp.loc.fp,datfile,sep=" "))
}else{
####Load datlist containing SNP,EXP (preloaded into matrixeqtl), and anno data 
  load(datfile)
}

###Function for cross validation

mat.train <- function(snpdat,expdat,train.indices,MEQTL.params){
  snpdat$ColumnSubsample(train.indices)
  expdat$ColumnSubsample(train.indices)
  with(MEQTL.params
    Matrix_eQTL_main(
      snps=snpdat,
      gene=expdat,
      cvrt=cvrt,
      output_file_name=paste0(output.file.name.tra,train.indices[1]),
      output_file_name.cis=paste0(output.file.name.cis,train.indices[1]),
      useModel=useModel,
      errorCovariance=errorCovariance,
      verbose=verbose,
      pvOutputThreshold=pvOutputThreshold.tra,
      pvOutputThreshold.cis=pvOutputThreshold.cis,
      snpspos=snpspos,
      genepos=genepos,
      cisDist=cisDist,
      pvalue.hist=pvalue.hist
    )
  )
  
}

samples <- datfile[["snps"]]$nCol()

train.indices <- chunk(samples,chunk.size=)

train.indices <- chunk(rep(1:513,9),n.chunks=9)


mapply(FUN=function(x,y)x[-y],)



test.indices <- chunk(1:513,chunk.size=57)







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


