args <- commandArgs(T)

require("reshape2")
require("tidyr")
require("dplyr")
require("Matrix")

raw_data <- read.table(args[1],header = F)
#raw_data <- read.table('K562_6sG_4sU_TFEA_cell_umi_gene.txt',header = F)
raw_data$gene <- gsub('--N','',raw_data$V3)
raw_data$gene <- gsub('--A','',raw_data$gene)
raw_data$gene <- gsub('--C','',raw_data$gene)

cells.keep <- raw_data %>% dplyr::distinct(V1,V2,gene) %>% group_by(V1) %>% dplyr::summarize(count=n()) %>% arrange(desccharacter 
inds <- as.numeric(args[2])
#inds <- 6000
if (length(cells.keep) > inds) {
  cells.keep2 <- head(cells.keep,inds)
}

raw_data <- raw_data %>% filter(V1 %in% cells.keep2) %>% droplevels

raw_data$type.A <- ''
if(TRUE %in% grepl('--A',raw_data$V3)){
  raw_data[grep('--A',raw_data$V3),]$type.A <- 'A'
}
raw_data$type.C <- ''
if(TRUE %in% grepl('--C',raw_data$V3)){
  raw_data[grep('--C',raw_data$V3),]$type.C <- 'C'
}

raw_data$type <- paste0(raw_data$type.A, raw_data$type.C)
count.gene <- count(raw_data, raw_data$V1, raw_data$gene, raw_data$type)

count.gene$gene <- paste0(count.gene$`raw_data$gene`, '--', count.gene$`raw_data$type`)
count.gene$gene <- as.factor(count.gene$gene)
count.gene$`raw_data$V1` <- as.factor(count.gene$`raw_data$V1`)
data.sparse = sparseMatrix(as.integer(count.gene$gene), as.integer(count.gene$`raw_data$V1`), x = count.gene$n)
colnames(data.sparse) = levels(count.gene$`raw_data$V1`)
rownames(data.sparse) = levels(count.gene$gene)
ord <- sort(colSums(data.sparse),decreasing = T)
data.sparse <- data.sparse[,names(ord)]
saveRDS(data.sparse,file=args[3])