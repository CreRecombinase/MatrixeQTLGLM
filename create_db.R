#Code to create expression and snp DBs
library(RSQLite)
library(MatrixEQTL)
library(sqldf)
library(reshape2)
library(fastmatch)
#USAGE Ctype Trans|Cis

oargs <- commandArgs(trailingOnly=T)
ctype <- oargs[1]
TC <- oargs[2]
setwd("/scratch/nwk2/pancan")
eqtl.file <- paste0("output/unimputed_",ctype,"_",TC,".txt")
eqtl.dat <- read.csv.sql(eqtl.file,sep="\t",header=T,eol="\n")

dat.file <- paste0(ctype,".Rdata")

load(dat.file)

snps <- as.matrix(snp.exp$snps)

exp <- as.matrix(snp.exp$gene)
rm(snp.exp)
gc()

snps <- snps[rownames(snps) %in% eqtl.dat$SNP,]
exp <- exp[rownames(exp) %in% eqtl.dat$gene,]

asnpnum <- apply(snps,1,function(x)sum(sort(tabulate(x+1),decreasing=T)[-1]))
agenenum <- apply(exp,1,function(x)sum(x!=0))

eqtl.dat$snpnum <- asnpnum[ fmatch(eqtl.dat$SNP,names(asnpnum))]
eqtl.dat$genenum <- agenenum[ fmatch(eqtl.dat$gene,names(agenenum))]

write.table(eqtl.dat,file=eqtl.file,sep="\t",col.names=T,row.names=F,quote=F,append=F)

dbfile <- paste0(ctype,"_",TC,".db")
db <- dbConnect(drv=dbDriver("SQLite"),dbname=dbfile)
bsnp <- melt(snps)
colnames(bsnp) <- c("Snp","Sample","Value")
dbGetQuery(db,"pragma journal_mode=memory")
dbGetQuery(db,"pragma synchronous=0")
dbGetQuery(db,"pragma cache_size=250000")
dbWriteTable(db,"snps",bsnp,overwrite=T,append=F)
#dbGetQuery(db,"pragma temp_store=1")
#dbGetQuery(db,"pragma temp_store_directory='.'")
#dbGetQuery(db,"create index ss on snps(Snp,Sample")

rm(bsnp)
gc()

bexp <- melt(exp)
colnames(bexp) <- c("Gene","Sample","Value")
dbWriteTable(db,"gene",bexp,row.names=F,overwrite=T,append=F)
dbDisconnect(db)
system(paste0("~/glm_eqtl/MatrixeQTLGLM/import-index.sh -x -g -d ",dbfile))
system(paste0("~/glm_eqtl/MatrixeQTLGLM/import-index.sh -x -s -d ",dbfile))
















