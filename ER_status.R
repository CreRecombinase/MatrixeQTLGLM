#Code to type ER status
library(sqldf)
exp.data <- read.csv.sql("C:/Users/nknoblau/Documents/R_WS/MatrixeQTLGLM/brca_RNAseq_expression.txt",sep="\t",header=T,eol="\n")
rownames(exp.data)<- exp.data[,1]
exp.data <- exp.data[,-1]
exp.data <- data.matrix(exp.data)
colnames(exp.data) <- gsub("_","-",colnames(exp.data))



plot(sort(exp.data["ESR1",]))

hist(exp.data["ESR1",])


#Train logistic regression model??

g.erstat <- read.delim("C:/Users/nknoblau/Documents/R_WS/MatrixeQTLGLM/clinical_patient_brca.txt",stringsAsFactors=F,header=T,fill=T)
colnames(exp.data)[!colnames(exp.data)%in% g.erstat$bcr_patient_barcode]

g.erstat <- g.erstat[ g.erstat$breast_carcinoma_estrogen_receptor_status %in% c("Positive","Negative"),]


g.erstat$ERSTAT <- ifelse(g.erstat$breast_carcinoma_estrogen_receptor_status=="Positive",1,0)

g.erstat <- g.erstat[ g.erstat$bcr_patient_barcode %in% colnames(exp.data),]

boxplot(split(exp.data["ESR1",g.erstat$bcr_patient_barcode],g.erstat$ERSTAT))

er.positive.seq <- exp.data[,g.erstat$bcr_patient_barcode[g.erstat$ERSTAT==1]]
er.negative.seq <- exp.data[,g.erstat$bcr_patient_barcode[g.erstat$ERSTAT==0]]


write.table(er.positive.seq,"brca_RNAseq_positive.txt",sep="\t",col.names=T,row.names=T,quote=F)
write.table(er.negative.seq,"brca_RNAseq_negative.txt",sep="\t",col.names=T,row.names=T,quote=F)

rnaseq.data <- read.csv.sql("brca_RNAseq_negative.txt",sep="\t",header=T,eol="\n")


