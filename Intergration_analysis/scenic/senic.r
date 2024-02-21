###This script computes transcription factors by using R package, 
###The first step of the SCENIC formal analysis is to calculate the correlation between transcription factors and each gene,
### which consumes a large amount of computational resources.Two strategies are recommended to deal with this step.When using GENIE3 to infer co-expression modules, 
###a small number of cells are randomly selected for calculation, and all cells are used for calculating the activity of regulons. The Python version of SCENIC is used to infer co-expression modules.


library(Seurat)
library(tidyverse)
library(patchwork)
library(SCENIC)
args<-commandArgs(T)
sample_name<-args[1]
type=args[2]
scenicOptions<-readRDS("scenicOptions.rds")
exprMat<-read.delim(paste0(sample_name,".txt"))
exprMat<-as.matrix(exprMat)
genesKept <- geneFiltering(exprMat, scenicOptions, 
              minCountsPerGene = 3 * 0.01 * ncol(exprMat), 
              minSamples = ncol(exprMat) * 0.01)
exprMat_filtered <- exprMat[genesKept, ]
saveRDS(exprMat_filtered,file="exprMat_filtered.rds")
exprMat_filtered<-as.matrix(exprMat_filtered)
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1)
runGenie3(exprMat_filtered_log, scenicOptions, nParts = 20)
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions,coexMethod=c("w005", "top50"))
exprMat_all <- as.matrix(exprMat_filtered)
exprMat_all <- log2(exprMat_all+1)
runSCENIC_3_scoreCells(scenicOptions, exprMat=exprMat_all)
runSCENIC_4_aucell_binarize(scenicOptions, exprMat=exprMat_all)







AUCmatrix <- readRDS("int/3.4_regulonAUC.Rds")
AUCmatrix <- AUCmatrix@assays@data@listData$AUC
AUCmatrix <- data.frame(t(AUCmatrix), check.names=F)
RegulonName_AUC <- colnames(AUCmatrix)
RegulonName_AUC <- gsub(' \\(','_',RegulonName_AUC)
RegulonName_AUC <- gsub('\\)','',RegulonName_AUC)
colnames(AUCmatrix) <- RegulonName_AUC
AUCmatrix[1:7,1:3]
celltype1<-subset(celltype1,orig.ident!="shang1-IV-2")
celltype1@meta.data=celltype1@meta.data[,c(1,2,3,4,5,6,7,19,20)]
head(celltype1)
x <- row.names(AUCmatrix) %>% gsub("\\.","-",.)%>%as.data.frame()
head(x)
row.names(AUCmatrix)<-x$.
AUCmatrix[1:10,1:3]
celltype1<-celltype1[,row.names(AUCmatrix)]
dim(celltype1)
scRNAauc <- AddMetaData(celltype1, AUCmatrix)
dim(scRNAauc@meta.data)
outfile1<-paste0(type,'_AUCmatrix.txt')
outfile2<-paste0(type,'_AUC.rds')

write.table(AUCmatrix,file=outfile1 , sep ="\t", row.names =TRUE, col.names =TRUE, quote =FALSE)
scRNAauc <- RunUMAP(scRNAauc, dims = 1:30, verbose = F)

scRNAauc <- FindNeighbors(scRNAauc, dims = 1:30, verbose = F)
scRNAauc <- FindClusters(scRNAauc, resolution = 0.5, verbose = F)
saveRDS(scRNAauc,file=outfile2)
