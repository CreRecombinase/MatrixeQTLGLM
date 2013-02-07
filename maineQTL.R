#First revision of Andy's pseudocode for unimputed BRCA SNP and Expression Data
#2/6/13
#NWK
library(sqldf)
library(plyr)
library(BatchExperiments)
library(MatrixEQTL)

snp.type <- "unimputed"
cancer.type <- "brca"
root.dir <- "/scratch/nwk2/mEQTL_ERpnc/glmEQTL/unimputed_brca/"
out.dir <- paste0(root.dir,"MEQTL_outputs/")

setwd(root.dir)
snp.filepath <- paste(snp.type,,cancer.type,"snp.txt",sep="_")
exp.filepath <- paste(cancer.type,"expression.txt",sep="_")

snp.loc.fp <- paste(snp.type,cancer.type,"snp_anno.txt",sep="_")
exp.loc.fp <- paste(cancer.type,"expression_anno.txt")

datfile <- "static.Rdata"

if(!file.exists(datfile)){
  shell(paste("/home/nwk2/glm_eqtl/MatrixeQTLGLM/load_static.R","F",snp.filepath,exp.filepath,snp.loc.fp,exp.loc.fp,datfile,sep=" "))
}else{
####Load datlist containing SNP,EXP (preloaded into matrixeqtl), and anno data 
  load(datfile)
}

###Function for cross validation

mat.train <- function(snpdat,expdat,train.indices,MEQTL.params){
  total.ids <- snpdat$nCol()
  kf <- sort(setdiff(train.indices,1:total.ids))
  snpdat$ColumnSubsample(train.indices)
  expdat$ColumnSubsample(train.indices)
  with(MEQTL.params
    Matrix_eQTL_main(
      snps=snpdat,
      gene=expdat,
      cvrt=cvrt,
      output_file_name=paste0(output.file.name.tra,kf,".txt"),
      output_file_name.cis=paste0(output.file.name.cis,kf,".txt"),
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

train.indices <- chunk(rep(1:samples,9),n.chunks=9)
#57 is a factor of 513, the number of samples
test.indices <- chunk(1:samples,chunk.size=57)


train.indices <- mapply(FUN=function(x,y)x[-y],train.indices,test.indices,SIMPLIFY=F)

MEQTL.params <- list(
  output.file.name.tra=paste0(out.dir,snp.type,"_",cancer.type,"_trans"),
  output.file.name.cis=paste0(out.dir,snp.type,"_",cancer.type,"_cis"),
  useModel=modelLINEAR,
  errorCovariance=SlicedData$new(),
  verbose=T,
  pvOutputThreshold.tra=1e-8,
  pvOutputThreshold.cis=1e-8,
  snpspos = datlist[["snp.anno"]],
  genepos = datlist[["exp.anno"]],
  cisDist=1e6
  pvalue.hist=F
)

m.dir <- tempfile(paste0("meqtl.res",cancer.type,"_",snp_type),tmpdir=out.dir)

MEQTL.reg <- makeRegistry(paste0("meqtl_reg_",cancer.type),file.dir=m.dir,packages="MatrixEQTL")

batchMap(MEQTL.reg,mat.train,train.indices=train.indices,more.args=list(snpdat=datlist[["snps"]],expdata=datlist[["gene"]],MEQTL.params=MEQTL.params))

submitJobs(MEQTL.reg)














s1<-sample(c(1:n),pts,replace=T)
for(i in 1:n){
  exp.train<-exp[s1!=i,]
  snp.train<-snp[s1!=i,]
  exp.test<-exp[s1==i,]M
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

