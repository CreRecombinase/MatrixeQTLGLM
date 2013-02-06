#!/home/nwk2/R-2.15.2/bin/Rscript
##Script to read in static data and save it to a Rdata file 
#nifty, as it reads from the command line
library(MatrixEQTL)
library(sqldf)

datlist <- list()

#Usage  <Use_MatrixEQTL_FileReader(T|F)> <SNP_File> <Expression_File > <SNP_Annotation_File> <Expression_Annotation_File> <Output_Rdata_Path>
#Returns a list in that order

#Helper Functions
load.data.matrix <- function(filepath){
  #Reads in data matrices
  rawdat <- read.csv.sql(filepath,sep="\t",header=T,eol="\n")
  rownames(rawdat)<- rawdat[,1]
  rawdat <- rawdat[,-1]
  rawdat <- data.matrix(rawdat)
  colnames(rawdat)<- gsub("_","-",colnames(rawdat))
  return(rawdat)
}
load.anno <- function(filepath){
  #Reads in Annotations
  read.csv.sql(filepath,sep="\t",header=T,eol="\n")
}

args <- commandArgs(trailingOnly=TRUE)

datlist$snp.anno <- load.anno(args[4])
datlist$exp.anno <- load.anno(args[5])


if(!all(file.exists(args[2:5]))){
  badfiles <- args[2:5][!file.exists(args[2:5])]
  stop(paste0("file not found: ",badfiles))
}
  


if(args[1]=="T"){
  datlist[["snps"]] <- SlicedData$new()
  datlist[["snps"]]$fileDelimiter <- "\t"
  datlist[["snps"]]$fileOmitCharacters <- "null"
  datlist[["snps"]]$fileSkipRows <- 1
  datlist[["snps"]]$fileSkipColumns <- 1
  datlist[["snps"]]$LoadFile(args[2])
  
  datlist[["gene"]] <- SlicedData$new(load.data.matrix(args[3])) 
}else{
  datlist[["snps"]] <- SlicedData$new(load.data.matrix(args[2]))
  datlist[["gene"]] <- SlicedData$new(load.data.matrix(args[3]))
}
print(ls())
save(datlist,file=args[6])
