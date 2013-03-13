library(plyr)
library("reshape2")


#Code for generating Gene Level Matrices of gene expression data from downloaded Data Matrices
cancer.type <- "COAD"
data.files <- "R:/Beck Lab/Level_2_SNP_Data/fc90ffc8-c37a-45fe-9427-075dcaf5c63a.tar/Expression-Genes/UNC__AgilentG4502A_07_3/Level_3/"
expfiles <- dir(data.files,full.names=T,pattern="*TCGA*")


allexp <- ldply(expfiles,.fun=read.table,header=T,sep="\t",colClasses=c("character","character","numeric"),quote="",nrows=17814,na.strings=c("null","NULL"))


nallexp <- dcast(allexp,...~barcode)
rownames(nallexp)<-nallexp[,1]
nallexp<-nallexp[,-1]

nbcodes <- substr(colnames(nallexp),start=1,stop=12)
nbcodes[1]

TN <- substr(colnames(nallexp),start=14,stop=16)

tnexp <- nallexp[,TN=="01A"]
colnames(tnexp) <- nbcodes[TN=="01A"]


rownames(tnexp) <- gsub("'","",rownames(tnexp))
outfile <- paste0("R:/Beck Lab/Level_2_SNP_Data/Level_3_Expression_Data/",cancer.type,"_Level_3_Expression.txt")
write.table(tnexp,file=outfile, sep="\t",col.names=NA,row.names=T,quote=F)


