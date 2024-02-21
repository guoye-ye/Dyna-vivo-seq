args <- commandArgs(T)
# args[1] <- 'K562_4sU_star_gene_exon_tagged_TagIntronic_clean_1000.TagTC.corrected_gene_cell_UMI_read.txt'
txt <- read.table(args[1], sep = '\t', na.strings = 'NA')
gene_mutation <- txt$V1
gene_mutation <- strsplit(as.character(gene_mutation), split = '--')
gene_mutation <- do.call(rbind,gene_mutation)
data <- cbind(gene_mutation,txt[,2:4])
tmp<-strsplit(args[1], '/', fixed = F)
len1<-length(tmp[[1]])
prefix<-tmp[[1]][len1]
prefix <- strsplit(prefix, '\\.', fixed = F)[[1]][1]
print(prefix)
filename <- paste0(prefix, '_tidy.txt')
write.table(data, file = filename, quote = F, sep = ' ', row.names = F, col.names = F)