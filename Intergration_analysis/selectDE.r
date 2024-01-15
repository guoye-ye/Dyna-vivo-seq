
library(dplyr)
library(Seurat)
library(tidyverse)


seurat_obj <- subset(seurat_obj,orig.ident!="IRI10min-1" & orig.ident!="IRI10min-2")
meta=seurat_obj@meta.data
meta$group=meta$orig.ident
meta$group=as.character(meta$group)
meta$orig.ident=as.character(meta$orig.ident)
meta <- meta %>% 
  mutate(group = case_when(
    group == "IRI30min-1" ~ "IRI30min",
    group == "IRI30min-2" ~ "IRI30min",
    group == "IRI50min-TSO1" ~ "IRI50min",
    group == "J-GZ-15h-1" ~ "J-GZ-1.5h",
    group == "J-GZ-15h-2" ~ "J-GZ-1.5h",
    group == "J-GZ-3h-1" ~ "J-GZ-3h",
    group == "J-GZ-3h-2" ~ "J-GZ-3h",
    group == "J-GZ-NC-1" ~ "J-GZ-NC",
    TRUE ~ group # 如果都不符合就返回原始值
  )
)
seurat_obj@meta.data=meta
split.list <- SplitObject(seurat_obj , split.by = "group")


all_norm_0_matrix<-as.data.frame(split.list[["J-GZ-NC"]]@assays[["RNA"]]@data)[row.names(matrix1),]
all_norm_1_matrix<-as.data.frame(split.list[["IRI30min"]]@assays[["RNA"]]@data)[row.names(matrix1),]
all_norm_2_matrix<-as.data.frame(split.list[["IRI50min"]]@assays[["RNA"]]@data)[row.names(matrix1),]
all_norm_3_matrix<-as.data.frame(split.list[["J-GZ-1.5h"]]@assays[["RNA"]]@data)[row.names(matrix1),]
all_norm_4_matrix<-as.data.frame(split.list[["J-GZ-3h"]]@assays[["RNA"]]@data)[row.names(matrix1),]
q_0<-apply(all_norm_0_matrix,1,mean)
q_1<-apply(all_norm_1_matrix,1,mean)
q_2<-apply(all_norm_2_matrix,1,mean)
q_3<-apply(all_norm_3_matrix,1,mean)
q_4<-apply(all_norm_4_matrix,1,mean)
all_merge_expr_mean1<-data.frame(J_GZ_NC=as.vector(q_0),IRI30min=as.vector(q_1),IRI50min=as.vector(q_2),J_GZ_1.5h=as.vector(q_3),J_GZ_3h=as.vector(q_4),row.names=row.names(matrix1))
write.table(all_merge_expr_mean1, file = 'DEmeanfortime.txt',sep ="\t", row.names =TRUE, col.names =TRUE, quote =FALSE)


pheatmap(select,scale="row",border=FALSE,fontsize_col=10, cluster_cols=FALSE,fontsize=11,cellwidth=30,cellheight=8,filename = "Pt3highDE.pdf")
write.table(allmeanmatpos, file = 'allmeanmatpos.txt',sep ="\t", row.names =TRUE, col.names =TRUE, quote =FALSE)
