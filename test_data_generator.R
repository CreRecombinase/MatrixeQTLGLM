#script for generating test data for mat.train
library(MatrixEQTL)
#number of patients
N=50
#number of SNPs 
S <- 20
#number of genes
G <- 10

#Chromosomes
C <- 23
#Chromosome size
CS <- 1000000

#Gene size
GS <- 10000

#SNP data
SNPdata <- matrix(sample(0:2,N*S,replace=T),S,N)
norm.SNP <- (SNPdata-rowMeans(SNPdata))/apply(SNPdata,MARGIN=1,FUN=sd)
norm.SNP <- t(norm.SNP)

expdata <- matrix(runif(G*N),G,N)
norm.exp <- expdata-rowMeans(expdata)

bcor <- norm.exp%*%norm.SNP

tmat <-SlicedData$new(SNPdata)
tmat$RowStandardizeCentered()

tmat <- as.matrix(tmat)
mean(tmat[1,])


SNPpos <- sample(seq(CS),S,replace=F)
SNPchrm <- sample(seq(C),S,replace=T)
genechrm <- sample(seq(C),G,replace=T)
genestart <- sample(seq(CS),G,replace=F)
geneend <- genestart+GS
snp.name <- paste("SNP_",seq(S),sep="")
gene.name <- paste("gene_",seq(G),sep="")

rownames(SNPdata)<- snp.name
rownames(expdata)<- gene.name
snp.anno <- data.frame(snp.name,SNPchrm,SNPpos,stringsAsFactors=F)
gene.anno <- data.frame(gene.name,genechrm,genestart,geneend)

snp.anno <- snp.anno[order(SNPchrm,SNPpos),]

gene.anno <- gene.anno[order(genechrm,genestart),]



snp.exp <- list(snps=SlicedData$new(SNPdata),gene=SlicedData$new(expdata))
save(snps.exp,file="test_snps_exp.Rdata")

annolist <- list(snp.anno=snp.anno,exp.anno=gene.anno)
save(annolist,file="test_anno.Rdata")

