##usage:Rscript mutation_rate_cmd_v2.R <file1> <file2> ...

library(ggplot2)
args<-commandArgs(T)
fileNum<-length(args)
print(paste0('file num: ',fileNum))


############# this function used to count the mutation ratio
CountMutation<-function(mutation,data,prefix){
  t_total<-data$V2[grepl('total_T',data$V1,perl = T)]
  t_rate<-(data$V2[grepl('^T_',data$V1,perl = T)]/t_total)[1:3]*100
  a_total<-data$V2[grepl('total_A',data$V1,perl = T)]
  a_rate<-(data$V2[grepl('^A_',data$V1,perl = T)]/a_total)[2:4]*100
  g_total<-data$V2[grepl('total_G',data$V1,perl = T)]
  g_rate<-(data$V2[grepl('^G_',data$V1,perl = T)]/g_total)[-2]*100
  c_total<-data$V2[grepl('total_C',data$V1,perl = T)]
  c_rate<-(data$V2[grepl('^C_',data$V1,perl = T)]/c_total)[-3]*100
  total_rate<-c(t_rate,a_rate,g_rate,c_rate)
  type<-data$V1[-c(1,5,6,7,11,13,16,19)]
  type<-gsub('_to_','>',type)
  group<-rep(prefix,12)
  mutation<-data.frame(rate=total_rate,type=type, group=group)
  return(mutation)
}

big_data<-data.frame()
prefix_list<-list()
for (i in 1:fileNum){
  assign(paste0('data',i),read.table(args[i]))
  assign(paste0('mutation',i),NULL)
  sample<-strsplit(args[i],split='/')
  tmp_l<-length(sample[[1]])
  sample<-sample[[1]][tmp_l]
  sample<-strsplit(sample,split='_')
  sample<-sample[[1]][1]
  print(sample)
  prefix_list[i]<-sample
  df<-CountMutation(mutation=get(paste0('mutation',i)),data=get(paste0('data',i)),prefix=sample)
  big_data<-rbind(big_data,df)
}
write.table(big_data,file="barbig_data.txt" , sep ="\t", row.names =TRUE, col.names =TRUE, quote =FALSE)



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
outfile<-paste0(outfile,'_muta_rate.pdf')


pdf(outfile,width = 12,height = 6)
ggplot(data=big_data, aes(x=type, y=rate,fill=group)) +
  geom_bar(stat="identity",color="black", position=position_dodge(0.9))+
  #geom_text(aes(label=total_rate), vjust=-0.5, size=5,position = position_dodge(1.1))+
  theme_minimal()+
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim = c(0,2.5))+
  ylab('Mutation rate / %')+
  xlab('Mutation type')+
  labs(title='The rate of different type of mutation ',hjust=0)+
  theme_bw()+theme(panel.border = element_rect(fill = NA, color = "black", size = 1.5, linetype = "solid"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  theme(legend.position = "top")+
  theme(legend.text = element_text(face='plain',size=16), # legend content size
        legend.title = element_blank(),# legend size
        #axis.line = element_line(size=1),   # Coordinate axis
        #text = element_text(family="Times"),# text font
        axis.text.y = element_text(size = 16,face='plain',color='black'),
        axis.text.x = element_text(size = 16,face='plain',color='black',angle=0, hjust=0.5, vjust=0.5),       
        axis.title = element_text(size=20,face="plain"),  
        plot.title = element_text(hjust=0.5,size=22,face=2))  
dev.off()