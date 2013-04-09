#Analysis of glm_eqtl network
library(plyr)
library(RSQLite)
library(reshape2)
library(sqldf)
library(boot)
library(fastmatch)
library(RcppArmadillo)
library(BatchJobs)

setwd("/scratch/nwk2/mEQTL_ERpnc/glmEQTL/brca_RNAseq/positive/overall/")
eqtl.file <- "unimputed_positive_trans.txt"
dbfile <- "overall_positive_linear.db"
pos.neg <- "Positive"
out.dir <- "."
n.chunks <- 500

trans.eqtl <- read.csv.sql(eqtl.file,header=T,sep="\t",eol="\n")



ntrans.eqtl <- trans.eqtl[ trans.eqtl$snpnum>12,]
ntrans.eqtl <- ntrans.eqtl[ ntrans.eqtl$genenum>20,]




a.chunks <- lapply(chunk(1:nrow(ntrans.eqtl),n.chunks=n.chunks),function(x)as.matrix(ntrans.eqtl[x,c("SNP","gene")]))


boot.ts <- function(t.chunks,dbname){
  db <- dbConnect(drv=dbDriver("SQLite"),dbname=dbname)
  gene.sql <- "select * from gene where gene=:gen"
  gene.mat <- dbGetPreparedQuery(db,gene.sql,data.frame(gen=unique(t.chunks[,"gene"])))
  gene.mat <- acast(gene.mat,Gene~Sample,value.var="Value")
  
  snp.sql <- "select * from snps where SNP=:sn"
  snp.mat <- dbGetPreparedQuery(db,snp.sql,data.frame(sn=unique(t.chunks[,"SNP"])))
  snp.mat <- acast(snp.mat,SNP~Sample,value.var="Value")
  
  srn <- rownames(snp.mat)
  grn <- rownames(gene.mat)
  dbDisconnect(db)
  
  gen.boot <- function(data,indices){
    data <- data[indices,]
    tf <- fastLmPure(data[,"snpval",drop=F],data[,"geneval"])
    return(tf$coefficients[1]/tf$stderr[1])
  }
  
  gen.t.stats <- function(SNP,gene){
    data <- cbind(snpval=snp.mat[fmatch(SNP,srn),],geneval=gene.mat[fmatch(gene,grn),])
    nboot <- boot(data=data,statistic=gen.boot,R=500)
    
    ti <- tryCatch(boot.ci(nboot,type="perc"),error=function(e)e)
    if(inherits(ti,"error")){
      return(NULL)
    }
    return(boot.ci(nboot,type="perc")$percent[4])
  }
  
  

    
  system.time(nats <- mdply(.data=t.chunks,.fun=gen.t.stats,.progress="time"))  
  return(nats)
  
}
m.dir <- tempfile(paste("boot_res_",pos.neg,sep=""),tmpdir=out.dir)
registry.name <- paste("boot_res_",pos.neg,sep="")
boot.reg <- makeRegistry(registry.name,file.dir=m.dir,packages=c("RSQLite","fastmatch","RcppArmadillo","boot","reshape2","plyr"))

batchMap(boot.reg,t.chunks =a.chunks,fun=boot.ts,more.args=list(dbname=dbfile))


#submitJobs(MEQTL.reg,resources=list(queue=args$QUEUE,memory=args$MEMORY,time=args$TIME,threads=1))




