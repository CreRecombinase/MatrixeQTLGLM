#Code to run a single overall MatrixEQTL (without covariates)
#3/8/13
#NWK
library(MatrixEQTL)

snp.type <- "imputed"
cancer.type <- "brca"
root.dir <- "/scratch/nwk2/mEQTL_ERpnc/imputed_eqtl/"
out.dir <- paste(root.dir,"output")

setwd(root.dir)


annofile <- "imputed_anno.Rdata"
snp.expdata <- "imputedSNP.genes.Rdata"
  
MEQTL.params <- list(
  output.file.name.tra=paste(out.dir,snp.type,"_",cancer.type,"_trans",sep=""),
  output.file.name.cis=paste(out.dir,snp.type,"_",cancer.type,"_cis",sep=""),
  useModel=modelLINEAR,
  verbose=T,
  pvOutputThreshold.tra=1e-8,
  pvOutputThreshold.cis=1e-8,
  cisDist=1e6,
  pvalue.hist=F
)
  

load(snp.expdata)
load(annofile)
with(MEQTL.params,Matrix_eQTL_main(
      snps=snp.exp$snps,
      gene=snp.exp$gene,
      output_file_name=paste(output.file.name.tra,i,".txt",sep=""),
      output_file_name.cis=paste(output.file.name.cis,i,".txt",sep=""),
      useModel=useModel,
      verbose=verbose,
      pvOutputThreshold=pvOutputThreshold.tra,
      pvOutputThreshold.cis=pvOutputThreshold.cis,
      snpspos=annolist$snp.anno,
      genepos=annolist$exp.anno,
      cisDist=cisDist,
      pvalue.hist=pvalue.hist
    )
  )
  
}