###This script is used to calculate the Pearson'r between the two samples. 
###
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(ggpubr)
#The next step here can be omitted depending on the actual situation。
data_pro <- function(data){
    data <- readRDS(data)
    data <- as.matrix(data)
    data <- CreateSeuratObject(counts=data, min.cells = 3)
    print(dim(data))
    data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
    data[["percent.rbp"]] <- PercentageFeatureSet(data, pattern = "^RP[SL]")
    filter_pbmc_filter <- subset(data, subset = nFeature_RNA > 500 & nFeature_RNA <5000 & percent.mt < 30) # 
    print(dim(filter_pbmc_filter))
    filter_pbmc_norm <- NormalizeData(filter_pbmc_filter, normalization.method = "LogNormalize", scale.factor = 10000)
    filter_pbmc_norm <- FindVariableFeatures(filter_pbmc_norm, selection.method = "vst", nfeatures = 2000)
    all.genes <- rownames(filter_pbmc_norm)
    filter_pbmc_norm <- ScaleData(filter_pbmc_norm, features = all.genes)
    return(filter_pbmc_norm)
}

CalPearson <- function(x1, x2, norm = FALSE){
    g1 <- rownames(x1)
    g2 <- rownames(x2)
    gene <- intersect(x= g1, y = g2)
    print(paste0('x1 gene number: ', length(g1)))
    print(paste0('x2 gene number: ', length(g2)))
    print(paste0('common gene number: ', length(gene)))
    if (norm == TRUE){
        x1 <- x1@assays$RNA@data
        x2 <- x2@assays$RNA@data
    }else{
        x1 <- x1@assays$RNA@counts
        x2 <- x2@assays$RNA@counts
    } 
    x1 <- x1[gene, ]
    x2 <- x2[gene, ]
    x1_mean <- rowMeans(x1)
    x2_mean <- rowMeans(x2)
    df <- data.frame(x1 = x1_mean, x2 = x2_mean)
    return(df)
}

myplot <- function(indata, inx, iny){
  nms <- names(indata)
  x <- nms[inx]
  y <- nms[iny]
  regression <- paste0(x, " ~ ", y)
  dat.lm <- lm(as.formula(regression), data = indata)
  r <- sprintf("italic(r) == %.2f",sqrt(summary(dat.lm)$r.squared))
  labels <- data.frame(r=r,stringsAsFactors = FALSE)
  
  ggplot(indata,aes(x=!!ensym(x), y=!!ensym(y)))+geom_point() + 
    geom_smooth(method = lm) + 
    labs(x=paste0(x," log(TP10K + 1)"),y=paste0(y," log(TP10K + 1)")) +
    geom_text(data=labels,mapping=aes(x = 0.1,y= 5.0,label=r),parse = TRUE,inherit.aes = FALSE,size = 5)+
    theme_bw()+theme_classic()+ ###用来移除top and rigth border
    theme(panel.grid.minor = element_blank())+
    theme(panel.grid.major = element_blank())+
    theme(axis.line.x=element_line(linetype=1,color="black",size=0.7),
          axis.line.y=element_line(linetype=1,color="black",size=0.7))+
    theme(plot.title = element_text(size=16,hjust = 0.5,vjust = 0.7,face="bold"))+
    theme(legend.title = element_blank(),
          legend.key.size = unit(0.5,'cm'),
          legend.text = element_text(size = 14),
          legend.position = "top")+
    theme(axis.title.x =element_text(size=14, face="bold"), 
          axis.title.y=element_text(size=14, face="bold"))+
    theme(axis.text = element_text(size = 12, face="bold"),
          axis.ticks.length.x=unit(-0.1, "cm"),
          axis.text.x = element_text(hjust=0.5,size = 12))+#angle=4
    coord_cartesian(xlim = c(0,5),ylim = c(0,5))
    
}

###total

sample <- 'shangxianc_total'
d1 <- data_pro('shang1_total.rds')  ##Change the name based on sample data
d2 <- data_pro('xia3_total.rds')
df <- CalPearson(d1, d2, norm = TRUE)
cov_pearson <- round(cov(df[,1], df[,2], method = 'pearson'),2)
outfile <- paste0(paste(sample, 'pearson', as.character(cov_pearson), sep='_'), '.pdf')

setwd("/R/R_code/yk")
pdf(outfile,width = 8, height = 8)
myplot(indata=df,inx=1,iny=2)
dev.off()




