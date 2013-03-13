#compare prediction of glmnet
library(sqldf)
library(MatrixEQTL)
install.packages("neuralnet")
setwd("C:/Users/nknoblau/Documents/R_WS/MatrixeQTLGLM/")
predict.mat <- read.table("predictionMatrix.csv",sep="\t",header=T,row.names=1,check.names=F)

predict.mat <- data.matrix(predict.mat)

nna.predict.mat <- predict.mat[apply(predict.mat,1,function(x)all(!is.na(x))),]
predict.mat <- nna.predict.mat

load("test/snp.exp.Rdata")

act.exp <- as.matrix(snp.exp$gene)
act.snp <- act.snp[unique(d.snps[[1]]),]


pred.genes <- intersect(rownames(nna.predict.mat),rownames(act.exp))

act.exp <- act.exp[pred.genes,colnames(predict.mat)]
act.snp <- as.matrix(snp.exp$snps)







LQsnps <- dbGetQuery(db,"select * from eqtls where Gene='LQK1' order by pValue")

head(LQsnps)

LQmat <- data.matrix(act.snp[unique(LQsnps[[1]]),])

plot(LQmat["rs963328",],act.exp["LQK1",])

t.exp <- as.vector(act.exp["LQK1",!is.na(act.exp["LQK1",])])
t.snp <- LQmat[,!is.na(act.exp["LQK1",])]


cv1 <- cv.glmnet(t(t.snp),t.exp,alpha=0.95,nfolds=10)
plot(cv1)

all.cors <- apply(LQmat,1,function(x)cor(x,t.exp,use="na.or.complete"))

 
which(all.cors==max(all.cors))

plot(LQmat[46,],act.exp["LQK1",])







se.mat <- (predict.mat-act.exp)^2
mse <- sqrt(apply(se.mat,1,sum,na.rm=T))
n.na <- apply(se.mat,1,function(x)sum(!is.na(x)))
mse <- mse/n.na

which(mse==min(mse))


any.nas <- apply(se.mat,MARGIN=1,function(x)!any(is.na(x)))

g.se.mat <- se.mat[any.nas,]
g.mse <- sqrt(apply(g.se.mat,1,sum))
which(g.mse==min(g.mse))




plot(predict.mat["POLR2K",],act.exp["POLR2K",])


cors <- numeric(nrow(predict.mat))
names(cors)<- rownames(predict.mat)
p.vals <- numeric(nrow(predict.mat))
names(p.vals) <- rownames(predict.mat)
rsquares <- numeric(nrow(predict.mat))
names(rsquares) <- rownames(predict.mat)
coeffs <- numeric(nrow(predict.mat))
names(coeffs)<- rownames(predict.mat)


for(i in 1:length(cors)){
  cors[i]<- cor(predict.mat[i,],act.exp[i,],use="na.or.complete",method="spearman")
  lmo <- lm(predict.mat[i,]~act.exp[i,],na.action=na.omit)
  f <- summary(lmo)$fstatistic
  rsquares[i]<- summary(lmo)$r.squared
  p <- pf(f[1],f[2],f[3])
  attributes(p)<- NULL
  p.vals[i]<- p
  coeffs[i]<- lmo$coefficients[2]
}

cors <- cors[rev(order(cors,na.last=F))]


coeffs <- coeffs[rev(order(coeffs,na.last=F))]
head(coeffs)

rsquares <- rsquares[rev(order(rsquares,na.last=F))]
head(rsquares)




