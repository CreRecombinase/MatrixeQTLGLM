#Rscript to melt data.frames from MatrixEQTL snp.exp lists
library(reshape2)
library(sqldf)
library(BBmisc)
library(RSQLite)

###USAGE  FILE <SNPS|GENE> (if SNPS) <CHUNK.SIZE> <DBFILE>

args <- list()

oargs <- commandArgs(trailingOnly=TRUE)
print(oargs)
  
args$FILE <- oargs[1]
args$SNPSGENE <- oargs[2]
if(args$SNPSGENE=="SNPS"){
  args$SNPCHUNKS=as.integer(oargs[3])
  args$DBFILE=oargs[4]
}else if(args$SNPSGENE=="GENE"){
  args$DBFILE=oargs[3]
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
db <- dbConnect(drv=dbDriver("SQLite"),dbname=args$DBFILE)
dbGetQuery(db,"pragma main.page_size=4096")
dbGetQuery(db,"pragma main.cache_size=50000")
dbGetQuery(db,"pragma synchronous=0")
dbGetQuery(db,"pragma journal_mode=memory")
if(dbExistsTable(db,"snps")){
  dbRemoveTable(db,"snps")
}



if(args$SNPSGENE=="GENE"){
  gene_melt <- melt(load.data.matrix(args$FILE))
  
  dbWriteTable(db,"genes",gene_melt,row.names=F,append=F,overwrite=T)
  
}else{
  snp.cols <- scan(args$FILE,what="character",nlines=1,sep="\t")[-1]
  conn <- file(args$FILE,open="r")
  ta <- read.table(conn,header=F,nrows=1)
  tsnp <- read.table(conn,header=F,nrows=args$SNPCHUNKS,sep="\t",stringsAsFactors=F)
  rownames(tsnp) <- tsnp[,1]
  tsnp <- tsnp[,-1]
  colnames(tsnp) <- snp.cols
  tsnp <- data.matrix(tsnp)
  ntsnp <- melt(tsnp)
  colnames(ntsnp) <- c("SNP","Sample","Value")
  dbWriteTable(db,"snps",ntsnp,row.names=F,append=T,overwrite=F)
  while(nrow(tsnp)==args$SNPCHUNKS){
    tsnp <- read.table(conn,header=F,nrows=args$SNPCHUNKS,sep="\t",stringsAsFactors=F)
    rownames(tsnp) <- tsnp[,1]
    tsnp <- tsnp[,-1]
    colnames(tsnp) <- snp.cols
    tsnp <- data.matrix(tsnp)
    ntsnp <- melt(tsnp)
    colnames(ntsnp) <- c("SNP","Sample","Value")
    dbWriteTable(db,"snps",ntsnp,row.names=F,append=T,overwrite=F)
  }
}
dbDisconnect(db)
