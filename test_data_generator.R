#script for generating test data for mat.train
library(MatrixEQTL)

test.data.directory <- "C:/Users/nknoblau/Documents/R_WS/MatrixeQTLGLM/testdata/"
setwd(test.data.directory)
#number of patients
N=50
#number of SNPs 
S <- 200
#number of genes
G <- 100

#number of transeQTLs
te <- 150
#number of ciseQTLs
ce <- 50

#Chromosomes
C <- 10
#Chromosome size
CS <- 1000

#Gene size
GS <- 100

#SNP data
SNPdata <- matrix(sample(0:2,N*S,replace=T),S,N)

expdata <- matrix(runif(G*N,min=-10,max=10),G,N)

SNPpos <- sample(seq(CS),S,replace=F)
SNPchrm <- sample(seq(C),S,replace=T)
genechrm <- sample(seq(C),G,replace=T)
genestart <- sample(seq(CS),G,replace=F)
geneend <- genestart+GS
snp.name <- paste("SNP_",seq(S),sep="")
gene.name <- paste("gene_",seq(G),sep="")


rownames(SNPdata)<- snp.name
rownames(expdata)<- gene.name
gene.eqtl <- sample(rownames(expdata),size=te,replace=T)
snp.eqtl <-  sample(rownames(SNPdata),size=te,replace=F)


colnames(SNPdata)<- paste0("Sample-",1:N)
colnames(expdata)<- paste0("Sample-",1:N)

for(i in 1:te){
  tk <- kmeans(expdata[gene.eqtl[i],],centers=3)
  ns <- match(tk$cluster,rank(tk$centers))-1
  SNPdata[snp.eqtl[i],] <- ns
}





SNPdata <- SNPdata[,sample(1:N,N,replace=F)]
expdata <- expdata[,sample(1:N,N,replace=F)]

snp.anno <- data.frame(snp.name,SNPchrm,SNPpos,stringsAsFactors=F)
gene.anno <- data.frame(gene.name,genechrm,genestart,geneend)

head(gene.eqtl,te)
for(i in 1:ce){
  snp.anno$SNPchrm[ snp.anno$snp.name %in% snp.eqtl[i]] <- gene.anno$genechrm[gene.anno$gene.name %in% gene.eqtl[i]]
  snp.anno$SNPpos[ snp.anno$snp.name %in% snp.eqtl[i]] <- sample(gene.anno$genestart[ gene.anno$gene.name %in% gene.eqtl[i]]:gene.anno$geneend[ gene.anno$gene.name %in% gene.eqtl[i]],size=1)
}


write.table(SNPdata,"testSNP.txt",col.names=NA,row.names=T,sep="\t",quote=F)
write.table(expdata,"testexp.txt",col.names=NA,row.names=T,sep="\t",quote=F)
write.table(snp.anno,"testSNPanno.txt",col.names=T,row.names=F,sep="\t",quote=F)
write.table(gene.anno,"testgeneanno.txt",col.names=T,row.names=F,sep="\t",quote=F)




snp.exp <- list(snps=SlicedData$new(SNPdata),gene=SlicedData$new(expdata))
save(snps.exp,file="test_snps_exp.Rdata")

annolist <- list(snp.anno=snp.anno,exp.anno=gene.anno)
save(annolist,file="test_anno.Rdata")

