#Rscript to melt data.frames from MatrixEQTL snp.exp lists
library(reshape2)
library(RSQLite)
library(MatrixEQTL)
library(BBmisc)

###USAGE  FILE <SNPS|GENE> (if SNPS) <SNPCHUNKS> (IF SNPS) <ROWNUM> <MELTFILE>

args <- list()

oargs <- commandArgs(trailingOnly=TRUE)
print(oargs)
  
args$FILE <- oargs[1]
args$SNPSGENE <- oargs[2]
if(args$SNPSGENE=="SNPS"){
  args$SNPCHUNKS=as.integer(oargs[3])
  args$ROWNUM=oargs[4]
  args$MELTFILE=oargs[5]
}else if(args$SNPSGENE=="GENE"){
  args$MELTFILE=oargs[3]
}else{
  stop(args)
}

load.data.matrix <- function(filepath){
  #Reads in data matrices
  rawdat <- read.csv.sql(filepath,sep="\t",header=T,eol="\n")
  rownames(rawdat)<- rawdat[,1]
  rawdat <- rawdat[,-1]
  rawdat <- data.matrix(rawdat)
  colnames(rawdat)<- gsub("_","-",colnames(rawdat))
  return(data.matrix(rawdat))
}


if(args$SNPSGENE=="GENE"){
  gene_melt <- melt(load.data.matrix(args$FILE))
  colnames(gene_melt)<- c("Gene","Sample","Value")
  write.table(gene_melt,file=args$MELTFILE,sep="\t",col.names=T,row.names=F,quote=F)
  
}else{
  snp.cols <- scan(args$FILE,what="character",nlines=1,sep="\t")[-1]
  snp.chunks <- floor(seq(from=1,to=floor(args$ROWNUM-(args$ROWNUM/args$SNPCHUNKS)),length.out=args$SNPCHUNKS))
  for(i in 1:length(snp.chunks)){
    tsnp <- read.csv.sql(args$FILE,header=F,sep="\t",eol="\n",skip=snp.chunks[i],
      sql=paste0("select * from file ", ifelse(i==length(snp.chunks),"",paste0("limit ",snp.chunks[i+1]-snp.chunks[i]))))
    rownames(tsnp) <- tsnp[,1]
    tsnp <- tsnp[,-1]
    colnames(tsnp) <- snp.cols
    tsnp <- data.matrix(tsnp)
    melt(tsnp)
    
    colnames(tsnp)<- c("Snp","Sample","Value")
    write.table(tsnp,file=args$MELTFILE,sep="\t",col.names=ifelse(i==1,T,F),row.names=F,quote=F,append=T)
  }
}