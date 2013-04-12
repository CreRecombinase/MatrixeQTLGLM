#!/home/nwk2/R-2.15.2/bin/Rscript
##Script to read in static data and save it to two Rdata files
#nifty, as it reads from the command line
library(MatrixEQTL)
library(sqldf)

annolist <- list()
snp.exp <- list()


#Usage <SNPEXP|ANNO>  (If Using arg SNPEXP) <Use_MatrixEQTL_FileReader(T|F)> <SNP_File> <Expression_File> <Output_Rdata_Path>
#Returns a list in that order
args <- list()
oargs <- commandArgs(trailingOnly=TRUE)
print(oargs)
args$EXPANNO <- oargs[1]
if(args$EXPANNO=="ANNO"){
  args$SNP <- oargs[2]
  args$EXP <- oargs[3]
  args$dpath <- oargs[4]
}else if(args$EXPANNO=="SNPEXP"){
  args$MEQTL <- oargs[2]
  args$SNP <- oargs[3]
  args$EXP <- oargs[4]
  args$dpath <- oargs[5]
}else{
  stop(args)
}



#Helper Functions
load.data.matrix <- function(filepath){
  #Reads in data matrices
  rawdat <- read.csv.sql(filepath,sep="\t",header=T,eol="\n")
  rawdat <- rawdat[!duplicated(rdat[,1]),]
  rownames(rawdat)<- rawdat[,1]
  rawdat <- rawdat[,-1]
  rawdat <- data.matrix(rawdat)
  colnames(rawdat)<- gsub("_","-",colnames(rawdat))
  return(data.matrix(rawdat))
}
load.anno <- function(filepath){
  #Reads in Annotations
  expargs <- scan(filepath,what="character",nlines=1,sep="\n")
  expargs <- strsplit(expargs,split="\t")[[1]]
  expargs <- expargs[-1]
  tf <- read.csv.sql(filepath,sep="\t",header=T,eol="\n")
  colnames(tf)<-expargs
  return(tf)
}


if(args$EXPANNO=="ANNO"){
annolist$snp.anno <- load.anno(args$SNP)
annolist$exp.anno <- load.anno(args$EXP)
save(annolist,file=args$dpath)

}else {
  if(args$MEQTL=="T"){
    snp.exp[["snps"]] <- SlicedData$new()
    snp.exp[["snps"]]$fileDelimiter <- "\t"
    snp.exp[["snps"]]$fileOmitCharacters <- "-1"
    snp.exp[["snps"]]$fileSkipRows <- 1
    snp.exp[["snps"]]$fileSkipColumns <- 1
    snp.exp[["snps"]]$LoadFile(args$SNP)
    snp.exp[["gene"]] <- SlicedData$new(load.data.matrix(args$EXP)) 
    sargs <- scan(args$SNP,what="character",nlines=1,sep="\n")
    sargs <- strsplit(sargs,split="\t")[[1]]
    sargs <- sargs[-1]
    expargs <- scan(args$EXP,what="character",nlines=1,sep="\n")
    expargs <- strsplit(expargs,split="\t")[[1]]
    expargs <- expargs[-1]
    snp.exp[["snps"]]$ColumnSubsample(match(expargs,sargs,nomatch=0))
    
    
  }else{
    snp.exp[["snps"]] <- SlicedData$new(load.data.matrix(args$SNP))
    snp.exp[["gene"]] <- SlicedData$new(load.data.matrix(args$EXP))
    sargs <- scan(args$SNP,what="character",nlines=1,sep="\n")
    sargs <- strsplit(sargs,split="\t")[[1]]
    sargs <- sargs[-1]
    expargs <- scan(args$EXP,what="character",nlines=1,sep="\n")
    expargs <- strsplit(expargs,split="\t")[[1]]
    expargs <- expargs[-1]
    snp.exp[["snps"]]$ColumnSubsample(match(expargs,sargs,nomatch=0))
  }
  save(snp.exp,file=args$dpath)
}

