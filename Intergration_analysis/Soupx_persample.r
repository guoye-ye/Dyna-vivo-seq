
library(Matrix)
library(SoupX)
library(Seurat)
library(DropletUtils)
args<-commandArgs(T)
fileNum<-length(args)
print(paste0('file num: ',fileNum))


process.danpi <- function(x){
  HL1<-x
  HL1<-as.matrix(HL1)
  HL1 <- CreateSeuratObject(counts=HL1, project='KN',min.cells = 3)
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




#设成function
groupinfo<-function(x){
  all <- x
  all <- CreateSeuratObject(all)
  all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
  all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 3000)
  all.genes <- rownames(all)
  all <- ScaleData(all, features = all.genes)
  all <- RunPCA(all, features = VariableFeatures(all), npcs = 30, verbose = F)
  all <- FindNeighbors(all, dims = 1:30)
  all <- FindClusters(all, resolution = 0.5)
  all <- RunUMAP(all, dims = 1:30)
}

runsoupx<-function(tod,matx){
  toc<-t(t(tod)[row.names(matx),])
  sc = SoupChannel(tod,toc)
  sc = setClusters(sc, setNames(matx$seurat_clusters, rownames(matx)))
  soupProf = data.frame(row.names = rownames(toc),est=rowSums(toc)/sum(toc),counts=rowSums(toc))
  sc = setSoupProfile(sc,soupProf)
  sc = autoEstCont(sc,verbose=FALSE)
  #out = adjustCounts(sc0,clusters=NULL,method='soupOnly',roundToInt=FALSE,verbose=1,tol=1e-3,pCut=0.01)
  out = adjustCounts(sc)
  return(out)
}
#过滤掉基因后需要保证基因数一致
#格式转换
#批量读取

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
  tod<-as.matrix(get(paste0('data',i)))
  KE<- process.danpi(get(paste0('data',i)))
  toc1<-KE@assays[["RNA"]]@counts
  all<-groupinfo(toc1)
  matx <- all@meta.data
  out<-runsoupx(tod,matx)
  prefix_list<-unlist(prefix_list)
  outfile<-NULL
  for (i in 1:fileNum){
    outfile<-paste0(prefix_list[i],'_total.rds')
    
  }
  saveRDS(out,file = outfile)
}



#HL1<-tod#分开判断稍仔细过滤
#HL1 <- CreateSeuratObject(counts=prefix_list[i], project='prefix_list[i]',min.cells = 3)





