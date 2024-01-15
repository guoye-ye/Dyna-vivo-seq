#this script used to  extract cell barcode from *_TC_matrix.rds (gene labeled with --T or --C)
#usega: Rscript extra_cb_from_TC_type_rds.R <*_TC_matrix.rd> <num_core_barcode>
args <- commandArgs(T)
matrix<-readRDS(args[1])
matrix<-as.matrix(matrix)
cb<-colnames(matrix)
cb<-as.data.frame(cb)
colnames(cb)<-'cells.keep2'
num_barcode<-args[2]
tmp<-strsplit(args[1], '/', fixed = F)
len1<-length(tmp[[1]])
prefix<-tmp[[1]][len1]
prefix<-gsub('_matrix.rds','',prefix)
outfile_tmp<-paste(prefix,'cb',num_barcode,sep="_")
outfile<-paste(outfile_tmp,'txt',sep='.')
write.table(cb,file=outfile,quote=F)