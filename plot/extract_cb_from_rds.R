args <- commandArgs(T)
#args[1] <- 'K562_4sU_TC_matrix.rds'
rds <- readRDS(args[1])
expression_matrix <- as.matrix(rds)
cb <- as.data.frame(colnames(expression_matrix))
cb <- as.data.frame(gsub('\n', '', cb$`colnames(expression_matrix)`))
tmp<-strsplit(args[1], '/', fixed = F)
len1<-length(tmp[[1]])
prefix<-tmp[[1]][len1]
prefix <- strsplit(prefix, '.', fixed = T)[[1]][1]
filename <- paste0(prefix, '_cb.txt')
write.table(cb, filename, quote = F, sep = '\t', col.names = F, row.names = F)