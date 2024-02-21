###This script is used to integrate multiple samples after initial filtering

library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(pheatmap)
library(harmony)
library(dplyr) 
library(Seurat)
library(tidyverse)
args<-commandArgs(T)
fileNum<-length(args)
print(paste0('file num: ',fileNum))



beforeharmony<-function(data,prefix){
  group_id<-prefix
  test.seu<-merge(x=get(paste0('data',1)),y=c(get(paste0('data',2)),get(paste0('data',3)),get(paste0('data',4)),get(paste0('data',5)),get(paste0('data',6)),get(paste0('data',7)),get(paste0('data',8))), add.cell.ids = group_id, project = 'big')
  print(paste0('dim:',dim(test.seu)))
  test.seu[["percent.mt"]]  <- PercentageFeatureSet(test.seu, pattern = "^mt-")
  test.seu[["percent.rbp"]] <- PercentageFeatureSet(test.seu, pattern = "^Rp[sl]")
  VlnPlot(test.seu, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbp"), ncol = 4)
  test.seu <- subset(test.seu, subset = nFeature_RNA > 500 & nFeature_RNA <6000 & percent.mt < 30)
  print(paste0('after filter:',dim(test.seu)))
  test.seu <- NormalizeData(object = test.seu, normalization.method = "LogNormalize", scale.factor = 1e4)
  test.seu <- FindVariableFeatures(object = test.seu, selection.method = 'vst', nfeatures = 2500)
  test.seu <- ScaleData(test.seu, features = rownames(test.seu),verbose = F)
  test.seu <- RunPCA(test.seu,features = VariableFeatures(object = test.seu), verbose = FALSE)
  test.seu <- FindNeighbors(test.seu, dims = 1:30, verbose = F)
  test.seu <- FindClusters(test.seu, resolution = 1.2, verbose = F)
  test.seu <- RunUMAP(test.seu, dims = 1:30, verbose = F)
  test.seu <- RunTSNE(test.seu, reduction = "pca", dims = 1:30)
  return(test.seu)
}
 


 prefix_list<-list()
for (i in 1:fileNum){
  assign(paste0('data',i),readRDS(args[i]))
  get(paste0('data',i))
  sample<-strsplit(args[i],split='/')
  tmp_l<-length(sample[[1]])
  sample<-sample[[1]][tmp_l]
  sample<-strsplit(sample,split='_')
  sample<-sample[[1]][1]
  print(sample)
  prefix_list[i]<-sample
}
  out<-list()  
  for (i in 1:fileNum){out[i]=get(paste0('data',i))}
  test.seu <- beforeharmony(data=get(paste0('data',i)),prefix=as.character(get("prefix_list")[1:fileNum]))
  prefix_list<-unlist(prefix_list)
  outfile<-NULL
  for (i in 1:fileNum){
  if(i==1){
    outfile<-paste(prefix_list[i],outfile,sep = '')
  }
  else{
    outfile<-paste(prefix_list[i],outfile,sep = '_')
  }
}
outfile1<-paste0(outfile,'_umapgroup1.2.pdf')
outfile2<-paste0(outfile,'_UMAP.cluster1.2.pdf')
outfile3<-paste0(outfile,'_total.rds')
outfile4<-paste0(outfile,'_umapgroup1.pdf')
outfile5<-paste0(outfile,'_UMAP.cluster1.pdf')
outfile6<-paste0(outfile,'_markers1.2.csv')
outfile7<-paste0(outfile,'_top_10_markers1.2.csv')




  pdf(file=outfile1)
  DimPlot(test.seu, label = F,pt.size = 1,group.by = "orig.ident")
  dev.off()
  pdf(file="tsnegroup.1.2.pdf")
  DimPlot(test.seu, reduction = "tsne",label=F, pt.size = 0.8,group.by = 'orig.ident')
  dev.off()
  pdf(file=outfile2)
  DimPlot(test.seu, reduction = "umap",label=T)
  dev.off()
markers <- FindAllMarkers(test.seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10_markers<-as.data.frame(markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC))
write.table(markers,'1.2_markers.txt',sep='\t',quote=FALSE,row.names=FALSE)
write.csv(top10_markers,file=outfile6)
write.csv(markers, file = outfile7)
saveRDS(test.seu,file=outfile3)

