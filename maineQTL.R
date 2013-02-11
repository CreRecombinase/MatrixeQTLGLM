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
out.dir <- paste(root.dir,"mEQTL_outputs/",sep="")

setwd(root.dir)
snp.filepath <- paste(snp.type,cancer.type,"snp.txt",sep="_")
exp.filepath <- paste(cancer.type,"expression.txt",sep="_")

snp.loc.fp <- paste(snp.type,cancer.type,"snp_anno.txt",sep="_")
exp.loc.fp <- paste(cancer.type,"expression_anno.txt")

annofile <- "annofile.Rdata"
snp.expdata <- "snp.exp.Rdata"
  
  
#if(!file.exists(datfile)){
#  shell(paste("/home/nwk2/glm_eqtl/MatrixeQTLGLM/load_static.R","F",snp.filepath,exp.filepath,snp.loc.fp,exp.loc.fp,datfile,sep=" "))
#}else{
####Load datlist containing SNP,EXP (preloaded into matrixeqtl), and anno data 
#load(annofile)
#}


###Function for cross validation
mat.train <- function(i,snp.exploc,anno.loc,train.indices,MEQTL.params){
  load(snp.exploc)
  load(anno.loc)
  total.ids <- snps.exp$snps$nCols()
  snps.exp$snps$ColumnSubsample(train.indices)
  snps.exp$gene$ColumnSubsample(train.indices)
  with(MEQTL.params,
    Matrix_eQTL_main(
      snps=snps.exp$snps,
      gene=snps.exp$gene,
      gene=snps.exp$exp,
      output_file_name=paste(output.file.name.tra,i,".txt",sep=""),
      output_file_name.cis=paste(output.file.name.cis,i,".txt",sep=""),
      useModel=useModel,
      verbose=verbose,
      pvOutputThreshold=pvOutputThreshold.tra,
      pvOutputThreshold.cis=pvOutputThreshold.cis,
      snpspos=annolist$snp.anno,
      genepos=annolist$exp.anno,
      cisDist=cisDist,
      pvalue.hist=pvalue.hist
    )
  )
  
}

col.command <- paste("head -1 ",exp.filepath," | awk '{print NF}'",sep="")

samples <- as.integer(system(col.command,intern=T))-1


train.indices <- chunk(rep(1:samples,9),n.chunks=9)
#57 is a factor of 513, the number of samples
test.indices <- chunk(1:samples,chunk.size=57)


train.indices <- mapply(FUN=function(x,y)x[-y],train.indices,test.indices,SIMPLIFY=F)

MEQTL.params <- list(
  output.file.name.tra=paste(out.dir,snp.type,"_",cancer.type,"_trans",sep=""),
  output.file.name.cis=paste(out.dir,snp.type,"_",cancer.type,"_cis",sep=""),
  useModel=modelLINEAR,
  verbose=T,
  pvOutputThreshold.tra=1e-8,
  pvOutputThreshold.cis=1e-8,
  cisDist=1e6,
  pvalue.hist=F
)

m.dir <- tempfile(paste("meqtl.res",cancer.type,"_",snp.type,sep=""),tmpdir=out.dir)
registry.name <- paste("meqtl_reg_",cancer.type,sep="")

MEQTL.reg <- makeRegistry(registry.name,file.dir=m.dir,packages="MatrixEQTL")

batchMap(MEQTL.reg,mat.train,train.indices=train.indices,i=1:length(train.indices),more.args=list(
  MEQTL.params=MEQTL.params,
  snp.exploc=snp.expdata,
  anno.loc=annofile))

submitJobs(MEQTL.reg)

