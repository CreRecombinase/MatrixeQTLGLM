#Code to generate list of lists to be fed into glm_predict.R
library(sqldf)
library(RSQLite)

args <- list()
##USAGE <RESULTS_DIR> <DATABSE_FILE> <KFOLD> <SAMPLE_NUM> <SAMPLE_FILE>

oargs <- commandArgs(trailingOnly=TRUE)
args$RESULTS_DIR <- oargs[1]
args$DATABASE_FILE <- oargs[2]
args$KFOLD <- oargs[3]
args$SAMPLE_NUM <- oargs[4]
args$SAMPLE_FILE <- oargs[5]



eqtl.files <- dir(args$RESULTS_DIR,pattern="*.txt",full.names=T)
kfolds <- as.integer(gsub(".+s([0-9]+).txt","\\1",eqtl.files))
cis.transs <- gsub(".+(cis|trans)[0-9]+.txt","\\1",eqtl.files)

write.db.file <- function(eqtl.file,cis.trans,kfold,dbfile){
  db <- dbConnect(drv=dbDriver("SQLite"),dbname=dbfile)
  eqtl <- read.csv.sql(eqtl.file,sep="\t",header=T,eol="\n")
  eqtl$Kfold<-kfold
  eqtl$cis.trans <- cis.trans
  dbWriteTable(db,name="eqtls",eqtl,row.names=F,append=T)
  dbDisconnect(db)
}

mapply(FUN=write.db.file,eqtl.file=eqtl.files,cis.trans=cis.transs,kfold=kfolds,MoreArgs=list(dbfile=args$DATABASE_FILE))






train.indices <- chunk(rep(1:args$SAMPLE_NUM,args$KFOLD),n.chunks=args$KFOLD)
test.indices <- chunk(1:args$SAMPLE_NUM,chunk.size=ceiling(args$SAMPLE_NUM/args$KFOLD))
train.indices <- mapply(FUN=function(x,y)x[-y],train.indices,test.indices,SIMPLIFY=F)


db <- dbConnect(drv=dbDriver("SQLite"),dbname=args$DATABSE_FILE)
dbSendQuery(db,"Create index sgk on eqtls(SNP,Gene,Kfold)")
dbSendQuery(db,"Create index gk on eqtls(Gene,Kfold)")

sample.cases <- scan(args$SAMPLE_CASES,what="character",sep="\n",nlines=1)
sample.cases <- strsplit(sample.cases,"\t")[[1]][-1]


test.df <- do.call("rbind",lapply(1:length(test.indices),function(x){
  data.frame(Sample=sample.cases[test.indices[[x]]],Kfold=x,Index=test.indices[[x]])
}))

train.df <-test.df <- do.call("rbind",lapply(1:length(train.indices),function(x){
  data.frame(Sample=sample.cases[train.indices[[x]]],Kfold=x,Index=train.indices[[x]])
}))

dbWriteTable(db,name="testSamples",test.df,row.names=F,overwrite=T,append=F)
dbWriteTable(db,name="trainSamples",train.df,row.names=F,overwrite=T,append=F)
dbSendQuery(db,"Create index k on testSamples(Kfold)")
dbSendQuery(db,"Create index ktr on trainSamples(Kfold)")

dbDisconnect(db)



