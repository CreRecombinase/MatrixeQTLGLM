#Code to generate list of lists to be fed into glm_predict.R

library(MatrixEQTL)
library(plyr)
library(doParallel)
library(sqldf)
if(Sys.info()['sysname']=="Windows"){
  snp.exploc <- "C:/Users/nknoblau/Documents/R_WS/MatrixeQTLGLM/test/snp.exp.Rdata"
  registerDoParallel(5)
  eqtl.base <- "C:/Users/nknoblau/Documents/R_WS/MatrixeQTLGLM/test2/unimputed_brca_trans"
  
  
  
}else{
  snp.exploc <- "/scratch/nwk2/mEQTL_ERpnc/glmEQTL/unimputed_brca/snp.exp.Rdata"
  registerDoParallel(12)
  eqtl.base <- "/scratch/nwk2/mEQTL_ERpnc/glmEQTL/unimputed_brca/57-fold/unimputed_brca_trans"
  save.dir <- "/scratch/nwk2/mEQTL_ERpnc/glmEQTL/unimputed_brca/test.train.snp.gene.Rdata"
  
}
eqtl.files <- paste0(eqtl.base,1:57,".txt")
sample.num <- 513
kfold <- 57



train.indices <- chunk(rep(1:sample.num,kfold),n.chunks=kfold)
test.indices <- chunk(1:sample.num,chunk.size=sample.num/kfold)
train.indices <- mapply(FUN=function(x,y)x[-y],train.indices,test.indices,SIMPLIFY=F)



aeqtls <- ldply(eqtl.files,function(x){
  tqtl <- read.csv.sql(x,header=T,eol="\n",sep="\t")
  knum <- gsub(".+trans([0-9]+).txt","\\1",x)
  tqtl$knum <- knum
  return(tqtl)
}
)

fsnp.genes <- split(aeqtls,aeqtls$knum)
snp.genes <- unlist(lapply(fsnp.genes,function(x)split(x$SNP,x$gene)),recursive=F)
rm(aeqtls,fsnp.genes)

load(snp.exploc)
snp.exp$snps <- t(as.matrix(snp.exp$snps))
snp.exp$gene <- t(as.matrix(snp.exp$gene))
snp.exp$gene <- snp.exp$gene[rownames(snp.exp$snps),]


knum <- as.integer(gsub("([0-9]+)\\..+","\\1",names(snp.genes)))



checkFun <- function(x){
  return(!all(apply(x,MARGIN=2,function(x)sum(sort(tabulate(x),decreasing=T)[-1]))<=2))
}


test.train.snp.gene <- llply(.data=1:length(snp.genes),function(i,snp.exp,snp.genes,train.indices,test.indices){
  gname <- gsub("[0-9]+.(.+)","\\1",names(snp.genes)[i])
  kn <- knum[i]
  train.snp <- snp.exp$snps[train.indices[[kn]],snp.genes[[i]],drop=F]
  train.exp <- snp.exp$gene[train.indices[[kn]],gname]
  badcols <- which(is.na(train.exp))
  if(length(badcols)>0){
    train.snp <- train.snp[-badcols,,drop=F]
    train.exp <- train.exp[-badcols]
  }
  test.snp <- snp.exp$snps[test.indices[[kn]],snp.genes[[i]],drop=F]
  if(checkFun(train.snp+1)){
    return(list(snp.train=train.snp,exp.train=train.exp,snp.test=test.snp,gn=gname))
  }
},snp.exp=snp.exp,snp.genes=snp.genes,train.indices=train.indices,test.indices=test.indices,.parallel=F,.paropts=list(multicombine=T,inorder=F),.progress="text")

