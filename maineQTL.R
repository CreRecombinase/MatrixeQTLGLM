#First revision of Andy's pseudocode for unimputed BRCA SNP and Expression Data
#2/6/13
#NWK
library(sqldf)
library(plyr)
library(BatchExperiments)
library(MatrixEQTL)

###USAGE maineQTL.R <out.files> <root.dir> <out-dir> <annofile> <snp.expfile> <samples> <fold-validation>
oargs <- commandArgs(trailingOnly=TRUE)
args <- list()
args$OUT.FILES <- oargs[1]
args$ROOT.DIR <- oargs[2]
args$OUT.DIR <- oargs[3]
args$ANNOFILE <- oargs[4]
args$SNP.EXPFILE <- oargs[5]
args$SAMPLES <- as.integer(oargs[6])
args$FOLD.VALIDATION <- as.integer(oargs[7])


root.dir <- args$ROOT.DIR
out.dir <- args$OUT.DIR

setwd(root.dir)


annofile <- args$ANNOFILE
snp.expdata <- args$SNP.EXPFILE
  


###Function for cross validation
mat.train <- function(i,snp.exploc,anno.loc,train.indices,MEQTL.params){
  load(snp.exploc)
  load(anno.loc)
  total.ids <- snp.exp$snps$nCols()
  snp.exp$snps$ColumnSubsample(train.indices)
  snp.exp$gene$ColumnSubsample(match(colnames(snp.exp$snps),colnames(snp.exp$gene)))
  with(MEQTL.params,
    Matrix_eQTL_main(
      snps=snp.exp$snps,
      gene=snp.exp$gene,
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

samples <- args$SAMPLES


train.indices <- chunk(rep(1:samples,args$FOLD.VALIDATION),n.chunks=args$FOLD.VALIDATION)


test.indices <- chunk(1:samples,chunk.size=ceiling(samples/args$FOLD.VALIDATION))


train.indices <- mapply(FUN=function(x,y)x[-y],train.indices,test.indices,SIMPLIFY=F)

MEQTL.params <- list(
  output.file.name.tra=paste(out.dir,args$OUT.FILES,"_trans",sep=""),
  output.file.name.cis=paste(out.dir,args$OUT.FILES,"_cis",sep=""),
  useModel=modelLINEAR,
  verbose=T,
  pvOutputThreshold.tra=1e-8,
  pvOutputThreshold.cis=1e-8,
  cisDist=1e6,
  pvalue.hist=F
)

m.dir <- tempfile(paste("meqtl.res",args$OUT.FILES,sep=""),tmpdir=out.dir)
registry.name <- paste("meqtl_reg_",args$OUT.FILES,sep="")

MEQTL.reg <- makeRegistry(registry.name,file.dir=m.dir,packages="MatrixEQTL")

batchMap(MEQTL.reg,mat.train,train.indices=train.indices,i=1:length(train.indices),more.args=list(
  MEQTL.params=MEQTL.params,
  snp.exploc=snp.expdata,
  anno.loc=annofile))



submitJobs(MEQTL.reg,resources=list(queue="short",memory=35000,time="3:00"))
Sys.sleep(35)
                      

