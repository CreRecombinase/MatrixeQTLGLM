#SNP Munging for GDAC SNP data
library(sqldf)
Hyb.name <- "C:/Users/nknoblau/Documents/R_WS/MatrixEQTL_ERpcn/SNP_names.txt"
SNP_anno <- "D:/Downloads/GenomeWideSNP_6.na32.annot.csv/nGWS.csv"

Hyb.dat <- read.csv.sql(Hyb.name,sep="\t",header=T,eol="\n")
anno.dat <- read.csv.sql(SNP_anno,sep=",",header=T,eol="\n")

Hyb.dat <- Hyb.dat[-1,]

anno.dat$Probe_Set_ID <- gsub("\"","", anno.dat$Probe_Set_ID)
anno.dat$dbSNP_RS_ID <- gsub("\"","", anno.dat$dbSNP_RS_ID)
anno.dat$Chromosome <- gsub("\"","", anno.dat$Chromosome)
anno.dat$Physical_Position <- gsub("\"","",anno.dat$Physical_Position)

rownames(anno.dat) <- anno.dat$Probe_Set_ID

n.anno <- anno.dat[Hyb.dat,"dbSNP_RS_ID"]
n.anno <- n.anno[!duplicated(n.anno)]

write(n.anno,"D:/Downloads/GenomeWideSNP_6.na32.annot.csv/n_anno.txt",sep="\n")

anno.dat <- anno.dat[ anno.dat$dbSNP_RS_ID %in% n.anno,]
anno.dat <- anno.dat[!duplicated(anno.dat$dbSNP_RS_ID),]

anno.dat$dbSNP_RS_ID[! anno.dat$dbSNP_RS_ID %in% n.anno]

anno.dat <- anno.dat[,c("dbSNP_RS_ID","Chromosome","Physical_Position")]
write.table(anno.dat,"D:/Downloads/GenomeWideSNP_6.na32.annot.csv/nGWS.txt",sep="\t",col.names=T,row.names=F,eol="\n",quote=F)

anno.dat <- anno.dat[n.anno,]

sum(duplicated(n.anno))

n.anno[duplicated(n.anno)]

RNAseq.dat <- read.csv.sql("C:/Users/nknoblau/Documents/R_WS/filter_bc_tumor_rnaseq.txt",sep="\t",header=T,eol="\n")

head(RNAseq.dat$Hybridization_REF,40)

RNAseq.dat$Hybridization_REF <- gsub("(.+)\\|.+","\\1",RNAseq.dat$Hybridization_REF)

RNAseq.dat <- RNAseq.dat[ RNAseq.dat$Hybridization_REF!="?",]
RNAseq.dat <- RNAseq.dat[!duplicated(RNAseq.dat$Hybridization_REF),]
colnames(RNAseq.dat) <- gsub("_","-",colnames(RNAseq.dat))


head(RNAseq.dat[,1:10])
write.table(RNAseq.dat,"C:/Users/nknoblau/Documents/R_WS/filter_bc_tumor_rnaseq.txt",sep="\t",col.names=T,row.names=F,quote=F)


RNA.anno <- read.csv.sql("C:/Users/nknoblau/Documents/R_WS/RNA_anno.txt",sep="\t",header=T,eol="\n")



head(RNA.anno)






