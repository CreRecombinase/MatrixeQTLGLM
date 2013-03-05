#Rscript to melt data.frames from MatrixEQTL snp.exp lists
library(reshape2)
library(RSQLite)
library(MatrixEQTL)
library(BBmisc)

###USAGE  <Rdata file> <SNPS|GENE> (if SNPS) <SNPCHUNKS> <MELTFILE>



oargs <- commandArgs(trailingOnly=TRUE)
print(oargs)
args <- list()
args$RDATA <- oargs[1]
args$SNPSGENE <- oargs[2]
if(args$SNPSGENE=="SNPS"){
  args$SNPCHUNKS=oargs[3]
  args$MELTFILE=oargs[4]
}else if(args$SNPSGENE=="GENE"){
  args$MELTFILE=oargs[3]
}else{
  stop(args)
}


load(args$RDATA)

if(args$SNPSGENE=="GENE"){
  gene_melt <- melt(as.matrix(snp.exp$gene))
  colnames(gene_melt)<- c("Gene","Sample","Value")
  write.table(gene_melt,file=args$MELTFILE,sep="\t",col.names=T,row.names=F,quote=F)
  
}else{
  snp.chunks <- chunk(1:nrow(snp.exp$snps),n.chunks=as.integer(args$SNPCHUNKS))
  for(i in 1:length(snp.chunks)){
    tsnp <- melt(as.matrix(snp.exp$snps)[snp.chunks[[i]],])
    colnames(tsnp)<- c("Snp","Sample","Value")
    write.table(tsnp,file=args$MELTFILE,sep="\t",col.names=ifelse(i==1,T,F),row.names=F,quote=F,append=T)
  }
}