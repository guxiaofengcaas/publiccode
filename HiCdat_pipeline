# run HiCdatR pipeline with R 

source("~/require.R")
library(pheatmap)
source("/HiCdatR/R/HiCdat-package.R")
chromSizes<-read.delim("/ref/9311/Bowtie2Index/genome.fa.fai",h=F)
f.internal.get.relevant.chromosomes <- function() {out <- chromSizes$V1;return(out) }
f.internal.get.chrom.sizes <- function() {out <- as.list(chromSizes$V2);names(out)<-chromSizes$V1;return(out)}
binSize=100000
rDir="./"

########### prepare matrix list ###########
samples <- c("Normal_rep1","Normal_rep2","HS_rep1","HS_rep2")
get_denseMatrix<-function(x){
  dataMatrix <-scan(paste0(x,"_",binSize,"_iced_dense.matrix"))
  dataMatrix <- matrix(dataMatrix, ncol = sqrt(length(dataMatrix)), byrow = TRUE)
  return(dataMatrix)
}

if(file.exists("binMatList.rdata")){
    load("binMatList.rdata")
}else{
    binMatList<-lapply(samples,get_denseMatrix)
    names(binMatList)<-samples
    save(binMatList,file="binMatList.rdata")
}
decay_data<-lapply(names(binMatList),function(x){
  data<-f.distance.decay(binMatList[[x]],binSize=binSize,rDir="./",distance=10e6,outfile=x,ncores=20)
  return(data)
})
names(decay_data)<-samples
f.plot_decay(decay_data,binSize = binSize)

############## Pairwise comparisons of Hi-C maps of different datasets

f.compare.samples.cor.difference(binMatList,binSize, rDir)


############## Correlation  of  differences

ccc <-
  f.plot.cor.difference(
    binMatList[["Normal_rep1"]],
    binMatList[["HS_rep1"]],
    binSize = binSize,
    rDir = rDir,
    randomizeDiff = FALSE,
    filterZero = F,
    outfile = "Normal_rep1_HS_rep1"
  )
f.Correlation_differences(chrID = "chr1",
                          data = ccc,
                          filename = "Normal_rep1_vs_HS_rep1.chr1.Correlation_differences.pdf")

############## Comparison  of  chromosome  1  Hi-C  maps

f.Comparison_3pic(
  filename = "HS_rep1_vs_Normal_rep1.chr1.pdf",
  binSize = binSize,
  dataMatrix1 = binMatList[["HS_rep1"]],
  dataMatrix2 = binMatList[["Normal_rep1"]]
)

