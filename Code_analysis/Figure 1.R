datapath='data/'
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
  write.table(DE_result,paste0('data/DEresult/mRNA_',compare_treat[k],'-',compare_control[k],'_DE.txt'),sep='\t',quote = F)
  
}
RBP_sig_all=data.frame()

for (i in 1:length(compare)) {
  compare0=compare[i]
  tempOutput=read.csv(paste0('DEresult/mRNA_',compare0,'_DE.txt'),sep='\t',header=T)
  tempOutput_RBP=tempOutput[which(rownames(tempOutput) %in% RBP),]
  tempOutput_RBP_sig=tempOutput_RBP %>% filter(adj.P.Val <= 0.05)
  tempOutput_RBP_sig$Gene=rownames(tempOutput_RBP_sig)
  if(nrow(tempOutput_RBP_sig)>0){
    RBP_sig_all0=cbind(tempOutput_RBP_sig,compare0)
    RBP_sig_all=rbind(RBP_sig_all,RBP_sig_all0)
  }
}

write.table(RBP_sig_all,'data/mRNA_DE_union.txt',sep='\t',row.names=F,quote=F)

RBP_sig_all=read.table('data/mRNA_DE_union.txt',sep='\t',header=T)
RBP_sig_allgene=unique(RBP_sig_all$Gene)
RBP_sig_0=unique(RBP_sig_all$Gene)
sig_all=data.frame()

for (i in 1:length(compare)) {
  compare0=compare[i]
  tempOutput=read.csv(paste0('data/DEresult/mRNA_',compare0,'_DE.txt'),sep='\t',header=T)
  tempOutput_sig=tempOutput %>% filter(adj.P.Val <= 0.05)
  tempOutput_sig$Gene=rownames(tempOutput_sig)
  if(nrow(tempOutput_sig)>0){
    sig_all0=cbind(tempOutput_sig,compare0)
    sig_all=rbind(sig_all,sig_all0)
  }
}

write.table(sig_all,'data/mRNA_allDE_union.txt',sep='\t',row.names=F,quote=F)

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

expdata_RBP_sig=expdata[RBP_sig_0,]
expdata_RBP_sig=expdata_RBP_sig[,names(sampletype2)]
expdata_RBP_sig_log=log2(expdata_RBP_sig+1)
library(ComplexHeatmap)
annot_df <- as.data.frame(sampletype2)
colnames(annot_df)='Disease_development'
colsmy_control = list(Disease_development= c( "Critical"="#99000D","Severe"="#FB6A4A",
                                              "Mild"="#FCBBA1", "Asymptomatic"="#FFF5F0"))

gap_col=as.numeric(table(annot_df))

cellInfo2_n3=c(sum(gap_col[1]),sum(gap_col[1:2]),sum(gap_col[1:3]))

top_anno <- HeatmapAnnotation(df=annot_df,col=colsmy_control)

SARS3=read.csv('data/SARS_RBP_intersect.csv',header=T)
SARS3_RBP=unique(SARS3$Gene)

labels_row=rownames(expdata_RBP_sig_log)
index_at=which(labels_row %in% c(SARS3_RBP,'CORO1A','DHX29','RPS5'))

row_anno = rowAnnotation(foo = anno_mark(at = index_at,
                                         labels = labels_row[index_at]))
library(circlize)
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
pdf('mRNA_heatmap_2.pdf',width=10,height=12)
print(a)
dev.off()

deresult=paste0(workdir,'data/DEresult/')
compare=c('Mild-Asymptomatic',
          'Severe-Asymptomatic',
          'Critical-Asymptomatic',
          'Severe-Mild',
          'Critical-Mild',
          'Critical-Severe')
dat_all=data.frame()
for (i in 1:length(compare)) {
  compare0=compare[i]
  tempOutput=read.csv(paste0('data/DEresult/mRNA_',compare0,'_DE.txt'),sep='\t',header=T)
tempOutput_RBP=tempOutput[which(rownames(tempOutput) %in% RBP),]
tempOutput_other=tempOutput[-which(rownames(tempOutput) %in% RBP),]
library(dplyr)
tempOutput_RBP_sig=tempOutput_RBP %>% filter(adj.P.Val <= 0.05 )
RBP_sig=nrow(tempOutput_RBP_sig)

tempOutput_RBP_nosig=tempOutput_RBP %>% filter(adj.P.Val > 0.05 )
RBP_nosig=nrow(tempOutput_RBP_nosig)

tempOutput_other_sig=tempOutput_other %>% filter(adj.P.Val <= 0.05 )
other_sig=nrow(tempOutput_other_sig)
tempOutput_other_nosig=tempOutput_other %>% filter(adj.P.Val > 0.05 )
other_nosig=nrow(tempOutput_other_nosig)

datlist=rbind(c(RBP_sig,RBP_nosig),c(other_sig,other_nosig))
rownames(datlist)=c('RBP','non-RBP')
colnames(datlist)=c('DE','non-DE')
fisher_result=fisher.test(datlist)
fisher_result_p=fisher_result$p.value
fisher_result_OR=fisher_result$estimate
fisher_result_CI=fisher_result$conf.int

dati=c(RBP_sig,RBP_nosig,other_sig,other_nosig,fisher_result_p,fisher_result_OR,fisher_result_CI)
dat_all=rbind(dat_all,dati)
}

colnames(dat_all)=c('RBP_sig','RBP_nosig','other_sig','other_nosig','fisher_result_p','fisher_result_OR','fisher_result_CI_1','fisher_result_CI_2')
dat_all=cbind(compare,dat_all)
dat_all=as.data.frame(dat_all)
dat_all$compare=factor(dat_all$compare,levels = dat_all$compare)
write.csv(dat_all,'mRNA_RBPfisher_padjust.csv',quote=F,row.names=F)

dat_all=read.csv('mRNA_RBPfisher_padjust.csv',header=T)
dat_all$compare=factor(dat_all$compare,levels = dat_all$compare)
dat_all_bar=dat_all[,1:5]
dat_all_bar=reshape2::melt(dat_all_bar)
library(stringr)
dat_all_bar$Type=str_replace_all(dat_all_bar$variable,c(
  'RBP_sig'='RBP',
  'RBP_nosig'='RBP',
  'other_sig'='Other',
  'other_nosig'='Other'
))

#自己计算比例
dat_all_bar2=dat_all_bar %>% 
  group_by(compare,Type) %>% 
  mutate(ratio=value/sum(value)) %>%
  arrange(compare)
dat_all_bar2$x=str_replace_all(dat_all_bar2$compare,c(
  'Mild-Asymptomatic'='1',
  'Severe-Asymptomatic'='2',
  'Critical-Asymptomatic'='3',
  'Severe-Mild'='4',
  'Critical-Mild'='5',
  'Critical-Severe'='6'
))
dat_all_bar2$x=as.numeric(dat_all_bar2$x)
dat_all_bar2_1=dat_all_bar2 %>% filter(Type == 'RBP')
dat_all_bar2_2=dat_all_bar2 %>% filter(Type == 'Other')
compare=c('Mild-Asymptomatic',
          'Severe-Asymptomatic',
          'Critical-Asymptomatic',
          'Severe-Mild',
          'Critical-Mild',
          'Critical-Severe')
Disease_development= pal_startrek('uniform')(4)[c(3,2,4,1)]

library(ggplot2)
p_raio=ggplot() + 
  geom_bar(data=dat_all_bar2_2,
           aes(x=x,y=ratio,fill=variable),
           stat='identity',position = 'stack',width=0.4)+
  geom_bar(data=dat_all_bar2_1,
           aes(x=x+0.4+0.01,y=ratio,fill=variable),
           stat='identity',position = 'stack',width=0.4)+
  scale_x_continuous(breaks = rep(1.2:6.2,1),
                     labels = compare)+
  scale_fill_manual(values = c("#CDBE70","#CD6839","#5CACEE","#104E8B"))+
  theme_classic()+
  geom_text(aes(x=1.2,y=1.02),label="**",angle=270)+
  geom_text(aes(x=2.2,y=1.02),label="****",angle=270)+
  geom_text(aes(x=3.2,y=1.02),label="ns",angle=270)+
  geom_text(aes(x=4.2,y=1.02),label="****",angle=270)+
  geom_text(aes(x=5.2,y=1.02),label="ns",angle=270)+
  geom_text(aes(x=6.2,y=1.02),label="****",angle=270)+
  labs(title = 'mRNA')+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )+
  coord_flip()
print(p_raio)

library(ggplot2)
dat_all_box=dat_all[,c(1,6:9)]
dat_all_box$Sig='Notsig'
dat_all_box[which(dat_all_box$fisher_result_p < 0.05),'Sig']='Sig'
dat_all_box$compare=factor(dat_all_box$compare,levels = compare)
dat_all_box=dat_all_box[order(dat_all_box$compare),]

p_OR=ggplot(data=dat_all_box, aes(x=compare,color=Sig)) + 
  geom_errorbar(aes(ymin = round(fisher_result_CI_1,4), ymax=round(fisher_result_CI_2,4)),
                width=0.2,
                position=position_dodge(0.9), 
                alpha = 1,
                size=1) + 
  theme_bw()+
  geom_point(aes(x=compare, y=round(fisher_result_OR,4)),pch=19,position=position_dodge(0.9),size=5) +
  scale_color_manual(values = c('#999999','#D6604D') )+
  theme(
    axis.text.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+  ylab('OR(95%CI)')+xlab('')+
  geom_hline(yintercept=1, linetype="dotted")+
  coord_flip()
print(p_OR)

library(patchwork)
pdf('Figure1B_new.pdf',width=8,height=8)
p_raio+p_OR+plot_layout(guides='collect',widths = c(1,1))
dev.off()