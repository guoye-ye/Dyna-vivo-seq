#Attention:This script does not distinguish between new RNA and old RNA
#usage:Rscript Generate_T_C_matrix_NoType_new_old_all_4cmd_v2.R <input> <total_RNA_num_barcode_file> <num_core_barcode> <sample_name> <data_type>

args <- commandArgs(T)

require("reshape2")
require("tidyr")
require("dplyr")
require("Matrix")
#setwd('/Users/ss/Documents/研究生/生信/yk/AKI/raw')
my.count1 <- read.table(args[1],h=F)
#my.count1 <- read.table("H1_star_gene_exon_tagged_TagIntronic_clean.TagTC.corrected_gene_cell_UMI_read_new.txt",h=F)

my.count1$V1 <- as.character(my.count1$V1)
my.count1$gene <- gsub("--C","",my.count1$V1)
my.count1$gene <- gsub("--T","",my.count1$gene)
cells.keep <- my.count1 %>% dplyr::distinct(V2,V3,gene) %>% group_by(V2) %>% dplyr::summarize(count=n()) %>% arrange(desc(count)) %>% .$V2 %>% as.character 
#cells.keep <- ifelse(length(cells.keep) > args[2],cells.keep[1:args[2]],cells.keep)# cells.keep[1:args[2]]
print("before trim")
print(length(cells.keep))
#inds <- as.numeric(args[2])#根据实际情况修改
#if (length(cells.keep) > inds) {
#  cells.keep2 <- head(cells.keep,inds)
#}

cells.keep2<-read.table(args[2],header = TRUE)
cells.keep2<-cells.keep2$cells.keep2
print("post trim")
print(length(cells.keep2))
#print(length(cells.keep))
my.count1 <- my.count1 %>% filter(V2 %in% cells.keep2) %>% droplevels
#print(dim(my.count1))

# check whether there are 5 columns
if (ncol(my.count1) !=5) {
  stop("Error! Please verify the count data frame!\n");
}
my.count2 <- my.count1 %>% arrange(gene,V2,V3)
my.count3 <- my.count2 %>% group_by(gene,V2) %>% dplyr::summarize(count=n())
my.count3$V2 <- as.factor(my.count3$V2)
my.count3$gene <- as.factor(my.count3$gene)
##将gene和cell barcode信息转换为坐标信息
data.sparse = sparseMatrix(as.integer(my.count3$gene), as.integer(my.count3$V2), x = my.count3$count)
colnames(data.sparse) = levels(my.count3$V2)
rownames(data.sparse) = levels(my.count3$gene)
ord <- sort(colSums(data.sparse),decreasing = T)
data.sparse <- data.sparse[,names(ord)]
dim(data.sparse)
num_core_barcode<-args[3]
sample_name<-args[4]
data_type<-args[5] #new ,old or total
outfile <- paste0(sample_name, '_TC_matrix_NoType_',data_type,sep = '_',num_core_barcode,'.rds')
saveRDS(data.sparse,outfile)