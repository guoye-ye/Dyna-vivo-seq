###This script is used to remove doublets in samples without TC labeling after running Soupx.
###usage:Rscript doublets_remove_per_sample.R <file.rds> 
#install
#library(remotes)
#remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
#rm(list=ls())
#setwd('~/R/R_code/IAA/intergration/vivo_seq')
#/usr/local/bin/R
library(Seurat)
library(Matrix)
library(DoubletFinder)
library(dplyr)
args<-commandArgs(T)
fileNum<-length(args)
print(paste0('file num: ',fileNum))



process.danpi <- function(data,prefix){
  HL1<-data
  HL1<- CreateSeuratObject(counts=HL1, project=prefix,min.cells = 3)
  HL1[["percent.mt"]] <- PercentageFeatureSet(HL1, pattern = "^mt-")
  HL1[["percent.rbp"]] <- PercentageFeatureSet(HL1, pattern = "^Rp[sl]")
  VlnPlot(HL1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rbp"), ncol = 4)
  HL1 <- subset(HL1, subset = nFeature_RNA > 200 & nFeature_RNA <6500) #
  print(dim(HL1))
  VlnPlot(HL1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  HL1 <- NormalizeData(HL1, normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
  HL1 <- FindVariableFeatures(HL1, selection.method = "vst", nfeatures = 2000, verbose = F)
  all.genes <- rownames(HL1)
  HL1 <- ScaleData(HL1, features = all.genes, verbose = F)
  HL1 <- RunPCA(HL1, features = VariableFeatures(object = HL1), verbose = F)
  ElbowPlot(HL1, ndim=50)
  HL1 <- FindNeighbors(HL1, dims = 1:30, verbose = F)
  HL1 <- FindClusters(HL1, resolution = 0.5, verbose = F)
  HL1 <- RunUMAP(HL1, dims = 1:30, verbose = F)
}


#TC <- process.danpi(TC)
rundoubletfinder<-function(x){
  TC<-x
  sweep.res.list <- paramSweep_v3(TC, PCs = 1:30, sct = F)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character()%>% as.numeric()
  #x<-unlist()
  #DoubletRate1 =x[2]*8*1e-6
  DoubletRate =dim(TC)[2]*8*1e-6
  homotypic.prop <- modelHomotypic(TC@meta.data[["seurat_clusters"]])
  #x = ncol(TC)
  nExp_poi <- round(DoubletRate*ncol(TC))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  TC <-doubletFinder_v3(TC,PCs = 1:30, pN = 0.25, pK = pK_bcmvn,nExp = nExp_poi.adj,reuse.pANN = F, sct = F)
  colnames(TC@meta.data)[9]<-"DS"
  DimPlot(TC, reduction = "umap", group.by = "DS")
  doubietinfo <- TC@meta.data
  #TC <- subset(TC, subset = DF.classifications_0.25_0.2_279 != "Doublet") 
  TC <- subset(TC, subset = DS != "Doublet") 
  #语法不对TC <- subset(TC, seurat_clusters == "7") 
  #TC@meta.data[["orig.ident"]]<-"prefix"
  #TC@project.name<-"KT_4SU"
  return(TC)
}

prefix_list<-list()
for (i in 1:fileNum){
  assign(paste0('data',i),readRDS(args[i]))
  sample<-strsplit(args[i],split='/')
  tmp_l<-length(sample[[1]])
  sample<-sample[[1]][tmp_l]
  sample<-strsplit(sample,split='_')
  sample<-sample[[1]][1]
  print(sample)
  prefix_list[i]<-sample
  TC <- process.danpi(data=get(paste0('data',i)),prefix=as.character(get("prefix_list")[i]))
  out<-rundoubletfinder(TC)
  prefix_list<-unlist(prefix_list)
  outfile<-NULL#
  for (i in 1:fileNum){
    outfile<-paste0(prefix_list[i],'_filter.rds') 
  }
  saveRDS(out,file = outfile)
}









