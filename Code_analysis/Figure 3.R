workdir='D:/Document/Postgraduate_study/Subject/RBP/'
datapath=paste0(workdir,'data/')
library(openxlsx)
library(ggplot2)
library(ggsci)
library(dylpr)
library(ggalluvial)
expdata=read.csv(paste0(datapath,'RBP_Gene_network.txt'),sep='\t',header=T)
expdata1=read.csv(paste0(datapath,'SARS_RBP_network.txt'),sep='\t',header=T)


colnames(expdata1) = c("SARS","RBP")
colnames(expdata) = c("RBP","RLDB")
expdata1 <- expdata1[1:19,]
expdata <- expdata[which(expdata$RBP %in% expdata1$RBP),]
expdata2 <- merge(expdata,expdata1,by="RBP")

expdata3 <- huzuobiao[which(huzuobiao$InteractionA %in% expdata2$RBP),]
rownames(expdata3) = c(1:6475)
colnames(expdata3) = c("RBP","x")

expdata4 <- huzuobiao[which(huzuobiao$InteractionB %in% expdata2$RBP),]
rownames(expdata4) = c(1:8793)
colnames(expdata4) = c("x","RBP")
expdata4$x1 <- expdata4$x
expdata4 <- expdata4[,-1]
colnames(expdata4) = c("RBP","x")

expdata5 <- merge(expdata3,expdata4,all = TRUE) %>% distinct()
freq <- as.data.frame(table(expdata5$x))
t <- as.data.frame(freq[which(freq$Freq >= 13),])

yuzhi13 <- expdata5[which(expdata5$x %in% t$Var1),]
yuzhi13 <- merge(yuzhi13,expdata2,by = "RBP")
yuzhi13 <- yuzhi13[,-3]

yuzhi13$SARS=factor(yuzhi13$SARS,levels=unique(yuzhi13$SARS))
yuzhi13=yuzhi13[order(yuzhi13$SARS),]

pdf('Interaction_13.pdf',width=10,height = 20)

ggplot(yuzhi13,
       aes(axis1 = SARS, axis2 = RBP, axis3 = x,fill=SARS)) +
  geom_flow(aes(fill=SARS),curve_type = "cubic",width = .2) +
  geom_stratum(width = 1/6, reverse = T,alpha = .5) +
  geom_text(stat = "stratum",aes(label = after_stat(stratum)),size=3,reverse=T) +
  scale_x_continuous(breaks = 1:3, labels = c("SARS", "RBP", "HUMAN")) +
  scale_fill_simpsons()+
  theme_minimal()+
  theme(panel.background = element_blank(),line = element_blank(),axis.text.y = element_blank())

dev.off()

SARS=read.csv('data/SARS-CoV-2.txt',sep='\t',header=T)
SARS=SARS[,c(1,3)]
extrainfo=read.xlsx('data/41422_2021_581_MOESM8_ESM.xlsx',sheet=2)
extrainfo$Gene1='SARS-CoV-2'
extrainfo=extrainfo[,c(ncol(extrainfo),1)]
colnames(extrainfo)=c('Gene1','Gene.2')
extrainfo=extrainfo[-which(extrainfo$Gene.2 %in% SARS$Gene.2),]
SARS=rbind(SARS,extrainfo)
RBP_union_info=read.xlsx('data/uniprot-yourlist.xlsx')
RBP_union_info=RBP_union_info[,c(1,2)]
colnames(RBP_union_info)[1]='Gene'
SARS2=merge(SARS,RBP_union_info,by.x='Gene.2',by.y='Entry',sort =F)
RBP_sig_all=read.csv('data/mRNA_DE_union.txt',sep='\t',header=T)
SARS3=merge(SARS2,RBP_sig_all,by='Gene',sort =F)
write.csv(SARS3,'data/SARS_RBP_intersect.csv',quote=F,row.names=F)

SARS3=read.csv('data/SARS_RBP_intersect.csv',header=T)
SARS3_2=SARS3 %>% filter(Gene1 != 'SARS-CoV-2')
head(SARS3_2)
humanP_humanP=read.csv('data/humanP_humanP.txt',sep='\t',header=T)
humanP_humanP_1=humanP_humanP %>% filter(InteractionA %in% SARS3_2$Gene) %>% distinct()
dim(humanP_humanP_1)
colnames(humanP_humanP_1)=c('RBP','Gene')
humanP_humanP_2=humanP_humanP %>% filter(InteractionB %in% SARS3_2$Gene) %>% distinct()
dim(humanP_humanP_2)
colnames(humanP_humanP_2)=c('Gene','RBP')

humanP_humanP_2=humanP_humanP_2[,c(2,1)]

allkk=rbind(humanP_humanP_1,humanP_humanP_2) %>% distinct()
dim(allkk)
kk=table(allkk$Gene)
head(kk)
for(i in 1:19){
  kk2=kk[kk >=i]
  print(paste(i,length(kk2)))
}

RBPselect=unique(SARS3_2$Gene)
Geneselect=names(kk[kk>=13])
genedup=intersect(RBPselect,Geneselect)
selectedgene=c(RBPselect,Geneselect)
selectedgene_type=c(rep('RBP',length(RBPselect)),rep('Gene',length(Geneselect)))
names(selectedgene_type)=selectedgene
selectedgene_type[which(names(selectedgene_type) == genedup)] = 'Both'
selectedgene_type2=selectedgene_type[!duplicated(names(selectedgene_type))]
selectedgene_type2=factor(selectedgene_type2,levels = c('RBP','Both','Gene'))
selectedgene_type2=sort(selectedgene_type2)

expdata=read.csv(paste0(datapath,'mRNA_FPKM_FinalVersion_20210415.txt'),sep='\t',header=T)
library(dplyr)
colnames(expdata)=lapply(colnames(expdata),function(a){
  b=unlist(strsplit(a,'_'))[1]
  return(b)
}) %>% unlist()
sampletype=read.csv(paste0(datapath,'Patientinfo.txt'),sep='\t',header=T)
sampletype2=sampletype$Status
names(sampletype2)=sampletype$PatientID
sampletype2=sampletype2[which(names(sampletype2) %in% colnames(expdata))]
sampletype2=factor(sampletype2,levels = c('Asymptomatic','Mild','Severe','Critical'))
sampletype2=sampletype2[order(sampletype2)]

expdata_RBP_sig=expdata[names(selectedgene_type2),]
expdata_RBP_sig=expdata_RBP_sig[,names(sampletype2)]
expdata_RBP_sig_log=log2(expdata_RBP_sig+1)
library(ComplexHeatmap)
annot_df <- as.data.frame(sampletype2)
colnames(annot_df)='Disease_development'
annot_df2 <- as.data.frame(selectedgene_type2)
colnames(annot_df2)='GeneClass'
colsmy_control = list(Disease_development= c( "Critical"="#99000D","Severe"="#FB6A4A",
                                              "Mild"="#FCBBA1", "Asymptomatic"="#FFF5F0"),
                      GeneClass = c('RBP' = '#F06574','Gene' = '#22CDF9','Both'='#97C763')
)
library(RColorBrewer)
library(pheatmap)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
p_heatmap =pheatmap::pheatmap(expdata_RBP_sig_log,scale = "row",fontsize = 10,breaks=bk,
                              annotation_row=annot_df2,annotation_colors=colsmy_control,
                              annotation_col = annot_df,
                              main='',cluster_rows=FALSE,cluster_cols = FALSE,
                              show_colnames = F ,
                              color=colorRampPalette(rev(c("red","white","blue")))(length(bk)),
)

pdf('E:/B-RBP_target_heatmap.pdf',width=8,height=12)
print(p_heatmap)
dev.off()

png('E:/B-RBP_target_heatmap.png',width=1000,height=1500)
print(p_heatmap)
dev.off()

expdata2=reshape2::melt(as.matrix(expdata))
colnames(expdata2)=c('Gene','PatientID','Expression')
expdata2=merge(expdata2,sampletype,by='PatientID')
Disease_development= pal_startrek('uniform')(4)[c(3,2,4,1)]

for(i in 1:length(Geneselect)){
  Geneselect0=Geneselect[i]
  df2=expdata2 %>% filter(Gene == Geneselect0)
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
    labs(title = paste0(Geneselect0,'_mRNA'))+
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
  
  pdf(paste0('C-mRNA_boxplot_RBPtarget_',Geneselect0,'.pdf'),width=4,height=8)
  print(p)
  dev.off()
  try(dev.off())
}