####### USAGE:Rscripts library_noLabeling_cmd_V002.R <infile1> [<infile2>] ··· 
# make sure your execute directory is 'library'

library(RColorBrewer)
library(ggplot2)
args<-commandArgs(T)

########## 1.data prepare function
data_prepare<-function(file){
  len<-length(strsplit(file,'/')[[1]])
  infile<-strsplit(file,'/')[[1]][len]
  data<-read.table(file,header = T)
  prefix<-gsub('_out_gene_exon_tagged_[0-9]+.dge.summary.txt','',infile,fixed = FALSE)
  print(paste0('======>: STAR prepare ', prefix, ' file'))
  len2<-length(strsplit(infile,'_')[[1]])
  num_core_barcode<-strsplit(infile,'_')[[1]][len2]
  num_core_barcode<-strsplit(num_core_barcode,'\\.')[[1]][1]
  group<-rep(prefix,dim(data)[2])
  data_all<-cbind(data,group)
  return(list(dataframe = data_all, prefix = prefix, num_core_barcode = num_core_barcode))
}

########## 1.library plot function

plot_option<-function(){
  option<-theme_bw()+theme_classic()+ ###用来移除top and rigth border
    theme(panel.grid.minor = element_blank())+
    theme(panel.grid.major = element_blank())+
    theme(axis.line.x=element_line(linetype=1,color="black",size=1),
          axis.line.y=element_line(linetype=1,color="black",size=1))+
    theme(plot.title = element_text(size=24,hjust = 0.5,vjust = 0.7,face="bold"))+
    theme(legend.title = element_blank(),
          legend.key.size = unit(0.7,'cm'),
          legend.text = element_text(size = 18, face=2),
          legend.position = "right")+
    theme(axis.title.x =element_text(size=20, face="bold"), 
          axis.title.y=element_text(size=20, face="bold"))+
    theme(axis.text = element_text(size = 18, face="bold"),
          axis.ticks.length.x=unit(0.1, "cm"),
          axis.text.x = element_text(hjust=0.5))#angle=45,
  return(option)
}


# UMI_reads
UMI_plot<-function(data){
  ylim<-max(data$NUM_TRANSCRIPTS)+100
  p<-ggplot(data,aes(x=NUM_GENIC_READS,y=NUM_TRANSCRIPTS,color=group))+
    theme_bw() +theme(panel.border = element_rect(fill = NA, color = "black", size = 1.5, linetype = "solid"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
        theme(axis.title.x =element_text(size=20), 
        axis.title.y=element_text(size=20))+
        theme(legend.title = element_text(color="black", # 修改图例的标题
                              size=20, 
                              face="plain"),
        legend.text = element_text(color="black", # 设置图例标签文字
                           size = 18, 
                           face = "plain"))+
    geom_point()+
    coord_cartesian(xlim = c(0,50000),ylim = c(0,ylim))+
    xlab("Reads/cell") +
    ylab("UMI/cell")+
    labs(title = '') +
    theme(axis.text.x = element_text(size = 18,color = "black", # 颜色
                                 face = "plain", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                 vjust = 0.5, # 位置
                                 hjust = 0.5),
        axis.text.y = element_text(size = 18, 
                                    # family = "myFont", 
                                    color = "black",
                                    face = "plain", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
          axis.ticks=element_line(size = 1.5))+
    scale_color_brewer(palette="Set2")
    return(p)
}

# Genes/Reads
Gene_plot<-function(data){
  ylim<-max(data$NUM_GENES)+100
  p<-ggplot(data,aes(x=NUM_GENIC_READS,y=NUM_GENES,color=group))+
    theme_bw() +theme(panel.border = element_rect(fill = NA, color = "black", size = 1.5, linetype = "solid"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
        theme(axis.title.x =element_text(size=20), 
        axis.title.y=element_text(size=20))+
        theme(legend.title = element_text(color="black", # 修改图例的标题
                              size=20, 
                              face="plain"),
        legend.text = element_text(color="black", # 设置图例标签文字
                           size = 18, 
                           face = "plain"))+
    geom_point()+
    coord_cartesian(xlim = c(0,60000),ylim = c(0,ylim))+
    xlab("Reads/cell") +
    ylab("Genes/cell")+
    labs(title = '')+
    theme(axis.text.x = element_text(size = 18,color = "black", # 颜色
                                 face = "plain", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                 vjust = 0.5, # 位置
                                 hjust = 0.5),
        axis.text.y = element_text(size = 18, 
                                    # family = "myFont", 
                                    color = "black",
                                    face = "plain", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
          axis.ticks=element_line(size = 1.5))+
    scale_color_brewer(palette="Set2")
  return(p)
}

file_num<-length(args)
big_df<-NULL
big_pre<-''
big_num<-''
for (i in 1:file_num){
  file<-args[i]
  tmp<-data_prepare(file)
  sub_df<-tmp$dataframe
  prefix<-tmp$prefix
  num_core_barcode<-tmp$num_core_barcode
  if (i == 1 ){
    big_df<-sub_df
    big_pre<-prefix
    big_num<-num_core_barcode
  }else{
    big_df<-rbind(big_df,sub_df)
    big_pre<-paste0(big_pre,'_',prefix)
    big_num<-paste0(big_num,'_',num_core_barcode)
  }
}

if ((file_num) == 1){
  outfile1<-paste0(prefix,'_Abundance_UMIs_reads_',num_core_barcode,'.pdf')
  outfile2<-paste0(prefix,'_Abundance_Genes_reads_',num_core_barcode,'.pdf')
}else{
  outfile1<-paste0(big_pre,'Abundance_UMIs_reads_',big_num,'.pdf')
  outfile2<-paste0(big_pre,'Abundance_Genes_reads_',big_num,'.pdf')
}

pdf(file=outfile1,width = 10,height = 10)
UMI_plot(big_df)
dev.off()

pdf(file=outfile2,width = 10,height = 10)
Gene_plot(big_df)
dev.off()

