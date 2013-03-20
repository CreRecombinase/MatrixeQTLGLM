#Code to generate list of lists to be fed into glm_predict.R
library(sqldf)
library(plyr)
library(RSQLite)
library(BBmisc)
library(doParallel)
args <- list()
##USAGE <RESULTS_DIR> <DATABSE_FILE> <KFOLD> <SAMPLE_NUM> <SAMPLE_FILE> <threads>

oargs <- commandArgs(trailingOnly=TRUE)
args$RESULTS_DIR <- oargs[1]
args$DATABASE_FILE <- oargs[2]
args$KFOLD <- as.integer(oargs[3])
args$SAMPLE_NUM <- as.integer(oargs[4])
args$SAMPLE_FILE <- oargs[5]
args$THREADS <- as.integer(oargs[6])

registerDoParallel(args$THREADS-1)

eqtl.files <- dir(args$RESULTS_DIR,pattern="*.txt",full.names=T)
kfolds <- as.integer(gsub(".+s([0-9]+).txt","\\1",eqtl.files))
cis.transs <- gsub(".+(cis|trans)[0-9]+.txt","\\1",eqtl.files)

write.db.file <- function(eqtl.file,cis.trans,kfold,dbfile){
  db <- dbConnect(drv=dbDriver("SQLite"),dbname=dbfile)
  dbGetQuery(db,"pragma busy_timeout=20000")
  dbGetQuery(db,"pragma main.page_size=4096")
  dbGetQuery(db,"pragma main.cache_size=50000")
  dbGetQuery(db,"pragma synchronous=0")
  dbGetQuery(db,"pragma main.journal_mode=WAL")
  dbGetQuery(db,"pragma temp_store=1")
  dbGetQuery(db,"pragma temp_store_directory='.'")
  eqtl <- read.csv.sql(eqtl.file,sep="\t",header=T,eol="\n")
  eqtl$Kfold<-kfold
  eqtl$CisTrans <- cis.trans
  dbWriteTable(db,name="eqtls",eqtl,row.names=F,append=T)
  dbDisconnect(db)
}
qtlargs <- data.frame(eqtl.file=eqtl.files,cis.trans=cis.transs,kfold=kfolds,dbfile=args$DATABASE_FILE)
print(head(qtlargs))

m_ply(.data=qtlargs,.fun=write.db.file,.parallel=F,.paropts=list(.inorder=F,.export="write.db.file",.packages=c("RSQLite","sqldf")),.progress="text",.inform=T)



db <- dbConnect(drv=dbDriver("SQLite"),dbname=args$DATABASE_FILE)

dbSendQuery(db,"Create index sgkc on eqtls(SNP,gene,Kfold,CisTrans)")
dbSendQuery(db,"Create index gkc on eqtls(gene,Kfold,CisTrans)")
#########get just eqtls that are in every cross validation
#Subset on genes that have SOME prediction in every cross validation
dbSendQuery(db,"create table ggenes as SELECT gene FROM (SELECT gene,count(distinct Kfold) from eqtls group by gene having count(distinct Kfold)>9)")
distinct.genes <- dbGetQuery(db,"select distinct gene from ggenes")
distinct.snps <- dbGetQuery(db,"select distinct SNP from eqtls,ggenes where eqtls.gene=ggenes.gene")



write.table(distinct.snps,file="eqtlsnps.txt",sep="\t",col.names=T,row.names=F,quote=F)
write.table(distinct.genes,file="eqtlgenes.txt",sep="\t",col.names=T,row.names=F,quote=F)






train.indices <- chunk(rep(1:args$SAMPLE_NUM,args$KFOLD),n.chunks=args$KFOLD)
test.indices <- chunk(1:args$SAMPLE_NUM,chunk.size=ceiling(args$SAMPLE_NUM/args$KFOLD))
train.indices <- mapply(FUN=function(x,y)x[-y],train.indices,test.indices,SIMPLIFY=F)



sample.cases <- scan(args$SAMPLE_FILE,what="character",sep="\t",nlines=1)[-1]



test.df <- do.call("rbind",lapply(1:length(test.indices),function(x){
  data.frame(Sample=sample.cases[test.indices[[x]]],Kfold=x,Index=test.indices[[x]])
}))

train.df<- do.call("rbind",lapply(1:length(train.indices),function(x){
  data.frame(Sample=sample.cases[train.indices[[x]]],Kfold=x,Index=train.indices[[x]])
}))

dbWriteTable(db,name="testSamples",test.df,row.names=F,overwrite=T,append=F)
dbWriteTable(db,name="trainSamples",train.df,row.names=F,overwrite=T,append=F)
dbSendQuery(db,"Create index k on testSamples(Kfold)")
dbSendQuery(db,"Create index ktr on trainSamples(Kfold)")

dbDisconnect(db)



