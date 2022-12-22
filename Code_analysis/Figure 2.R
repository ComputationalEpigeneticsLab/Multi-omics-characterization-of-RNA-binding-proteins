expdata=read.csv(paste0(datapath,'Protein_FinalVersion_20210420_5pm.txt'),sep='\t',header=T)
rownames(expdata)=expdata[,1]
expdata=expdata[,-1]
sampletype=read.csv(paste0(datapath,'Patientinfo.txt'),sep='\t',header=T)
sampletype2=sampletype$Status
names(sampletype2)=sampletype$PatientID
library(dplyr)
colnames(expdata)=lapply(colnames(expdata),function(a){
  b=unlist(strsplit(a,'_'))[1]
  return(b)
}) %>% unlist()
sampletype2=sampletype2[which(names(sampletype2) %in% colnames(expdata))]
index_cleargene=which(apply(expdata,1,function(v){return((sum(v==0)/length(v))>=0.5)}))
if(length(index_cleargene)!=0){
  expdata2=expdata[-index_cleargene,]
}else{
  expdata2=expdata
}
expdata3=log2(expdata2 + 0.01)

compare_treat=c('Mild','Severe','Critical','Severe','Critical','Critical')
compare_control=c('Asymptomatic','Asymptomatic','Asymptomatic','Mild','Mild','Severe')

for(k in 1:length(compare_treat)){
  control_sample=names(sampletype2[sampletype2==compare_control[k]])
  treat_sample=names(sampletype2[sampletype2==compare_treat[k]])
  
  control<-rowMeans(expdata2[,control_sample])
  treat<-rowMeans(expdata2[,treat_sample])
  foldchange<- treat/control
  foldchange2<-as.data.frame(foldchange)
  foldchange2<-cbind(foldchange2,log2(foldchange))
  pvalue_all=c()
  for(i in 1:nrow(expdata3)){
    control=as.numeric(expdata3[i,control_sample])
    treat=as.numeric(expdata3[i,treat_sample])
    result=wilcox.test(control,treat)
    pvalue=result$p.value
    pvalue_all=c(pvalue_all,pvalue)
  }
  names(pvalue_all)=rownames(expdata3)
  padjust_all=p.adjust(pvalue_all,method = 'fdr')
  
  DE_result=cbind(foldchange2,pvalue_all,padjust_all)
  colnames(DE_result)=c('FC','logFC','P.Value','adj.P.Val')
  write.table(DE_result,paste0('data/DEresult/Protein_',compare_treat[k],'-',compare_control[k],'_DE.txt'),sep='\t',quote = F)
}

library(dplyr)
library(openxlsx)
RBPlist=read.xlsx(paste0(datapath,'RBP.xlsx'))
RBP=RBPlist$RBP.name
RBP=RBP[!is.na(RBP)]

RBP_sig_all=data.frame()

for (i in 1:length(compare)) {
  compare0=compare[i]
  tempOutput=read.csv(paste0('data/DEresult/Protein_',compare0,'_DE.txt'),sep='\t',header=T)
  tempOutput_RBP_sig=tempOutput %>% filter(adj.P.Val <= 0.05)
  tempOutput_RBP_sig$Gene=rownames(tempOutput_RBP_sig)
  if(nrow(tempOutput_RBP_sig)>0){
    RBP_sig_all0=cbind(tempOutput_RBP_sig,compare0)
    RBP_sig_all=rbind(RBP_sig_all,RBP_sig_all0)
  }
}

RBP_sig_0=unique(RBP_sig_all$Gene)
expdata=read.csv(paste0(datapath,'Protein_FinalVersion_20210420_5pm.txt'),sep='\t',header=T)
rownames(expdata)=expdata[,1]
expdata=expdata[,-1]
sampletype=read.csv(paste0(datapath,'Patientinfo.txt'),sep='\t',header=T)
sampletype2=sampletype$Status
names(sampletype2)=sampletype$PatientID

library(dplyr)
colnames(expdata)=lapply(colnames(expdata),function(a){
  b=unlist(strsplit(a,'_'))[1]
  return(b)
}) %>% unlist()
sampletype2=sampletype2[which(names(sampletype2) %in% colnames(expdata))]
sampletype2=factor(sampletype2,levels=c('Asymptomatic','Mild','Severe','Critical'))
sampletype2=sort(sampletype2)
expdata=expdata[,names(sampletype2)]

expdata_RBP_sig=expdata[RBP_sig_0,]
expdata_RBP_sig=expdata_RBP_sig[,names(sampletype2)]
expdata_RBP_sig_log=log2(expdata_RBP_sig+1)
compare=c('Mild-Asymptomatic',
          'Severe-Asymptomatic',
          'Critical-Asymptomatic',
          'Severe-Mild',
          'Critical-Mild',
          'Critical-Severe')
RBP_sig_Protein=read.table('data/Protein_allDE_union.txt',sep='\t',header=T)
RBP_sig_mRNA=read.table('data/mRNA_allDE_union.txt',sep='\t',header=T)

intersect_all=data.frame(Gene=character(),compare=character(),FC_status=character())
for(i in 1:length(compare)){
  compare00=compare[i]
  RBP_sig_Protein0=RBP_sig_Protein %>% filter(compare0 == compare00)
  RBP_sig_mRNA0=RBP_sig_mRNA %>% filter(compare0 == compare00)
  RBP_sig_Protein0_up = RBP_sig_Protein0[which(RBP_sig_Protein0$logFC > 0),]
  RBP_sig_mRNA0_up = RBP_sig_mRNA0[which(RBP_sig_mRNA0$logFC > 0),]
  k0_1=intersect(RBP_sig_Protein0_up$Gene,RBP_sig_mRNA0_up$Gene)
  if(length(k0_1) > 0){
    k1=cbind(k0_1,compare00,'Up')
  }else{k1=NA}
  RBP_sig_Protein0_down = RBP_sig_Protein0[which(RBP_sig_Protein0$logFC < 0),]
  RBP_sig_mRNA0_down = RBP_sig_mRNA0[which(RBP_sig_mRNA0$logFC < 0),]
  k0_2=intersect(RBP_sig_Protein0_down$Gene,RBP_sig_mRNA0_down$Gene)
  if(length(k0_2) > 0){
    k2=cbind(k0_2,compare00,'Down')
  }else{k2=NA}
  kk=rbind(k1,k2)
  colnames(kk)=c('Gene','compare','FC_status')
  kk=as.data.frame(na.omit(kk))
  intersect_all=rbind(intersect_all,kk)
}
library(ComplexHeatmap)
annot_df <- as.data.frame(sampletype2)
colnames(annot_df)='Disease_development'
colsmy_control = list(Disease_development= c( "Critical"="#99000D","Severe"="#FB6A4A",
                                              "Mild"="#FCBBA1", "Asymptomatic"="#FFF5F0")
)

top_anno <- HeatmapAnnotation(df=annot_df,col=colsmy_control)

labels_row=rownames(expdata_RBP_sig_log)
index_at=which(labels_row %in% c(intersect_all$Gene,'RPS5'))

row_anno = rowAnnotation(foo = anno_mark(at = index_at,
                                         labels = labels_row[index_at]))
library(circlize)
library(RColorBrewer)
col_fun2 = circlize::colorRamp2(c(-2, 0, 2), c("navy", "yellow", "firebrick3"))
col_fun2=colorRampPalette(rev(brewer.pal(11, "PiYG")))(11)
col_fun2=circlize::colorRamp2(c(-2, 0, 2),c('#FC8D59','#FFFFBF','#91BFDB'))
expdata_RBP_sig_log2=as.matrix(t(scale(t(expdata_RBP_sig_log))))
a<-Heatmap(expdata_RBP_sig_log2,
           cluster_rows = TRUE,
           cluster_columns = FALSE,
           show_column_names = FALSE,
           show_row_names = FALSE,
           top_annotation = top_anno,
           column_title = NULL,
           row_names_side = "left",
           row_names_gp = gpar(fontsize = 7),
           right_annotation = row_anno,
           heatmap_legend_param = list(title = "")
)
pdf('Protein_heatmap_2.pdf',width=10,height=18)
print(a)
dev.off()
pdf('Protein_heatmap.pdf',width= 7, height=3)
print(p)
dev.off()

dat_all=data.frame()
for (i in 1:length(compare)) {
  compare0=compare[i]
  tempOutput=read.csv(paste0('data/DEresult/mRNA_',compare0,'_DE.txt'),sep='\t',header=T)
  tempOutput_RBP=tempOutput[which(rownames(tempOutput) %in% RBP),]
  tempOutput_RBP$compare=compare0
  tempOutput_RBP$Symbol=rownames(tempOutput_RBP)
  tempOutput_RBP_sig=tempOutput_RBP %>% filter(P.Value < 0.05)
  dat_all=rbind(dat_all,tempOutput_RBP_sig)
}

write.csv(dat_all,'mRNA_RBP_sig_in6.csv',quote=F,row.names=F)
RBP_sig_num=as.data.frame(sort(table(dat_all$Symbol),decreasing=T))
write.csv(RBP_sig_num,'mRNA_RBP_sig_num.csv',quote=F,row.names=F)
dat_all=read.csv('mRNA_RBP_sig_in6.csv',header=T)
RBP_sig_num=read.csv('mRNA_RBP_sig_num.csv',header=T)
SARS3=read.csv('data/SARS_RBP_intersect.csv',header=T)
SARS3_RBP=unique(SARS3$Gene)

RBP_sig_top=RBP_sig_num %>% filter(Freq > 4)
RBP_sig_top=as.character(RBP_sig_top$Var1)

RBP_sig_top=c(RBP_sig_top,c('CORO1A','DHX29','RPS5'),SARS3_RBP) %>% unique()

datapath=paste0(workdir,'data/')
expdata=read.csv(paste0(datapath,'mRNA_FPKM_FinalVersion_20210415.txt'),sep='\t',header=T)
sampletype=read.csv(paste0(datapath,'Patientinfo.txt'),sep='\t',header=T)
sampletype2=sampletype$Status
names(sampletype2)=sampletype$PatientID

library(dplyr)
colnames(expdata)=lapply(colnames(expdata),function(a){
  b=unlist(strsplit(a,'_'))[1]
  return(b)
}) %>% unlist()
sampletype2=sampletype2[which(names(sampletype2) %in% colnames(expdata))]

RBPlist=read.xlsx(paste0(datapath,'RBP.xlsx'))
RBP=RBPlist$RBP.name
RBP=RBP[!is.na(RBP)]

expdata2=reshape2::melt(as.matrix(expdata))
colnames(expdata2)=c('Gene','PatientID','Expression')

expdata2=merge(expdata2,sampletype,by='PatientID')
Disease_development= pal_startrek('uniform')(4)[c(3,2,4,1)]
for(i in 1:length(RBP_sig_top)){
  RBP_sig_top0=RBP_sig_top[i]
  df2=expdata2 %>% filter(Gene == RBP_sig_top0)
  df2$Status=factor(df2$Status,levels=c('Asymptomatic','Mild','Severe','Critical'))
  df2=df2[order(df2$Status),]
  
  comparelist=list(c('Mild','Asymptomatic'),c('Severe','Asymptomatic'),c('Critical','Asymptomatic'),
                   c('Severe','Mild'),c('Critical','Mild'),
                   c('Critical','Severe'))
  library(ggpubr)
  library(ggsignif)
  p=ggplot(data=df2,aes(x=Status,y=log2(Expression),color=Status))+
    stat_boxplot(geom="errorbar",width=0.15,position = position_dodge(0.9))+
    geom_boxplot(position = position_dodge(0.9))+
    geom_jitter(aes(color=Status),size=1)+
    labs(title = paste0(RBP_sig_top0,'_mRNA'))+
    scale_color_manual(values = Disease_development)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=0,hjust = 1,colour = 'black',size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          legend.position = 'none',
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(color = 'black',
                                    size=10,
                                    hjust=0.5))+
    stat_compare_means(comparisons=comparelist,label = 'p.signif')
  
  pdf(paste0('data/RBP_exp_boxplot/mRNA/mRNA_boxplot_',RBP_sig_top0,'.pdf'),width=4,height=8)
  print(p)
  dev.off()
}

p=ggplot(data=df2,aes(x=Status,y=log2(Expression),color=Status))+
  stat_boxplot(geom="errorbar",width=0.15,position = position_dodge(0.9))+
  geom_boxplot(position = position_dodge(0.9))+
  geom_jitter(aes(color=Status),size=1)+
  labs(title = paste0(RBP_sig_top0,'_mRNA'))+
  scale_color_brewer(palette='Set1',direction =-1)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=0,hjust = 1,colour = 'black',size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        legend.position = 'none',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(color = 'black',
                                  size=10,
                                  hjust=0.5))+
  stat_compare_means(comparisons=comparelist,label = 'p.signif')