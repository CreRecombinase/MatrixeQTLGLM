#Comparison of RNAseq vs Gene Expression Array
library(sqldf)
library(limma)
RNA.seq.file <- "C:/Users/nknoblau/Documents/R_WS/filter_bc_tumor_rnaseq.txt"

rna.seq <- read.csv.sql(RNA.seq.file,sep="\t",header=T,eol="\n")

array.file <- "C:/Users/nknoblau/Documents/R_WS/MatrixEQTL_ERpcn/nLevel3.txt"
array.dat <- read.csv.sql(array.file,sep="\t",header=T,eol="\n")

rownames(rna.seq)<-rna.seq[,1]
rna.seq <- rna.seq[,-1]

rownames(array.dat) <- array.dat[,1]
array.dat <- array.dat[,-1]

inter.genes <- intersect(rownames(array.dat),rownames(rna.seq))
inter.cases <- intersect(colnames(array.dat),colnames(rna.seq))

array.dat <- array.dat[inter.genes,inter.cases]
rna.seq <- rna.seq[inter.genes,inter.cases]

array.dat <- data.matrix(array.dat)
rna.seq <- data.matrix(rna.seq)
rna.seq <- log(rna.seq+.Machine$double.eps)

norm.array <- array.dat-rowMeans(array.dat)

norm.seq <- rna.seq-rowMeans(rna.seq)







bad.rows <- apply(norm.array,1,function(x)all(is.na(x)))
norm.seq <- norm.seq[!bad.rows,]
norm.array <- norm.array[!bad.rows,]

coeffs <- numeric(nrow(norm.array))
p.vals <- numeric(nrow(norm.seq))
for(i in 1:nrow(norm.array)){
  tfit <- summary(lm(norm.seq[i,]~norm.array[i,]))
  coeffs[i] <- tfit$coefficients[2,1]
  f <- tfit$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  p.vals[i] <- p
  
}





