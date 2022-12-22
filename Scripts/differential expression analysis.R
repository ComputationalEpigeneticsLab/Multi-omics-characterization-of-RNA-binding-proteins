workdir='E:/Subject/RBP/'
setwd(workdir)
datapath=paste0(workdir,'data/')
deresult=paste0(workdir,'DEresult/')
library(dplyr)
library(openxlsx)
library(ggplot2)
library(ggsci)
library(openxlsx)
RBPlist=read.xlsx(paste0('data/RBP.xlsx'))
RBP=RBPlist$RBP.name
RBP=RBP[!is.na(RBP)]

humanP_humanP=read.csv('data/humanP_humanP.txt',sep='\t',header=T)

compare=c('Mild-Asymptomatic',
          'Severe-Asymptomatic',
          'Critical-Asymptomatic',
          'Severe-Mild',
          'Critical-Mild',
          'Critical-Severe')
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
mpg=as.data.frame(table(sampletype2))
mpg$sampletype2=factor(mpg$sampletype2,levels=c('Asymptomatic','Mild','Severe','Critical'))
mpg=mpg[order(mpg$sampletype2),]

Disease_development= pal_startrek('uniform')(4)[c(3,2,4,1)]
g <- ggplot(mpg, aes(sampletype2,Freq,fill=sampletype2))
g=g +  geom_col()+ theme_bw()+
  geom_text(aes(label=Freq), 
            position = position_dodge2(width = 0.8, preserve = 'single'), 
            vjust = -0.5, hjust = 0.5)+
  guides(fill=guide_legend(title='SampleType'))+scale_fill_manual(values=Disease_development)+
  xlab('Disease_development')+ylab('')+ggtitle('mRNA')
pdf('data/SampleType_count_mRNA.pdf')
print(g)
dev.off()

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

#wilcox.test
##treat_sample is more serious
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
  #false discovery rate
  padjust_all=p.adjust(pvalue_all,method = 'fdr')
  
  DE_result=cbind(foldchange2,pvalue_all,padjust_all)
  colnames(DE_result)=c('FC','logFC','P.Value','adj.P.Val')
  write.table(DE_result,paste0('data/DEresult/mRNA_',compare_treat[k],'-',compare_control[k],'_DE.txt'),sep='\t',quote = F)}

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