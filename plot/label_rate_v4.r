###This script was used for calculating the proportion of metabolically labeled newly synthesized single-cell RNA.It should be noted that this script can only be used for rough estimation of labeling rates 
###and is applicable for assessing metabolic labeling levels before data integration and filtering.
###usage:Rscript label_rate_cmd_v4.R <file1> <file2> ...
###file:*TC_mutation_count.txt
library(ggplot2)
args <- commandArgs(T)
fileNum<-length(args)
print(paste0('file num: ',fileNum))

fraction_calculate <- function(data){
  for(i in 1:dim(data)[1]){
    count_T <- data[i,"T"]
    count_C <- data[i,"C"]
    data[i,4] <- count_C + count_T
  }
  for(i in 1:dim(data)[1]){
    fraction <- data[i,"C"]/data[i,"V4"]
    data[i,5] <- round(fraction*100,4)
  }
  return(data)
}

big_data<-data.frame()
prefix_list<-list()
for (i in 1:fileNum){
  assign(paste0('data',i),read.table(args[i],header = T, sep = ','))
  sample<-strsplit(args[i],split='_')
  sample<-sample[[1]][1]
  print(sample)
  prefix_list[i]<-sample
  df<-fraction_calculate(data=get(paste0('data',i)))
  df[,6] <- sample
  df<-df[c("X","T","C","V4","V5","V6")]
  big_data <- rbind(big_data,df[,5:6])
}
colnames(big_data) <- c('fraction_of_labeled_transcripts_per_cell','type')


#### 如果要筛选总umi数大于多少就用以下代码
#T13$total<-T13$T+T13$C
#T14$total<-T14$T+T14$C
#T13<-subset(T13,total>300)
#T14<-subset(T14,total>300)

#T14<-subset(T14, V5>1)

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
outfile<-paste0(outfile,'_label_ratio_transcript.pdf')


pdf(outfile,width = 20,height = 10)
ggplot(big_data, aes(x = type, y = fraction_of_labeled_transcripts_per_cell)) +
  theme_bw()+
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
  geom_jitter(color="black",size=0.89,alpha=0.2,width=0.22,height=0.1)+
  geom_boxplot(width=0.3)+
  scale_fill_brewer(palette='Set1')+
  ylab('labeling ratio')+
  xlab('group')+
  theme(plot.title = element_text(size=18,hjust = 0.5,vjust = 0.5))+
  theme(axis.title.x =element_text(size=16), 
        axis.title.y=element_text(size=16))+
  theme(axis.text = element_text(size = 16))
  
dev.off()