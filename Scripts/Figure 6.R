
'''R
Rscript drug_optimization.R
'''


SARS3=read.csv('data/SARS_RBP_intersect.csv',header=T)
SARS3_RBP=unique(SARS3$Gene)
drug=read.csv('data/SARS_RBP_drugname.txt',sep='\t',header=T)#potential drug
drug=drug[order(drug$Z_score),]

#Drug-Target
drug_target=read.csv('data/drug_humanP.txt',sep='\t',header=T)
colnames(drug_target)[2]='Target'
drug2=merge(drug,drug_target,by.x='Drug',by.y='DrugID')
drug2=drug2[order(drug2$Z_score),]

humanP_humanP=read.csv('data/humanP_humanP.txt',sep='\t',header=T)
humanP_humanP_1=humanP_humanP %>% filter(InteractionA %in% drug2$Target)
humanP_humanP_2=humanP_humanP %>% filter(InteractionB %in% drug2$Target)
humanP_humanP_2=humanP_humanP_2[,c(2,1)]
humanP_humanP2=rbind(humanP_humanP_1,humanP_humanP_2)
colnames(humanP_humanP2)=c('Target','Symbol')

sarscov=read.csv('data/SARS-CoV-2.txt',sep='\t',header=T,check.names=F)
transid=read.csv('data/SARS-CoV-2_uniprotid_transsymbol.txt',sep='\t',header=T,check.names=F)
sarscov2=merge(sarscov,transid,by.x='Gene-2',by.y='From')
sarscov2=sarscov2[,c('Gene1','To')]
colnames(sarscov2)=c('Virus','Symbol')
extrainfo=read.xlsx('data/41422_2021_581_MOESM8_ESM.xlsx',sheet=2)
extrainfo$Gene1='SARS-CoV-2'
extrainfo=extrainfo[,c(ncol(extrainfo),2)]
colnames(extrainfo)=c('Virus','Symbol')
extrainfo=extrainfo[-which(extrainfo$Symbol %in% sarscov2$Symbol),]
sarscov_all=rbind(sarscov2,extrainfo)
head(drug2)
drug2_test=drug2[order(drug2$Z_score),]
length(unique(drug2_test$Name))

drug2_test=drug2[order(drug2$Z_score),]
drug2_test=drug2_test[1:60,]
drug2_test2=merge(drug2_test,humanP_humanP2,by='Target')
drug2_test2=drug2_test2 %>% filter(Symbol %in% SARS3_RBP)
drug2_test2=drug2_test2 %>% dplyr::select(Name,Target,Symbol)
colnames(drug2_test2)=c('Drug','Target','SARS_RBP')
data=drug2_test2
data$flow<-1:nrow(data)
data2<-reshape2::melt(data,id='flow')
variable<-summary(data2$variable)
data2$link<-1
library(randomcoloR)
mycolor<-distinctColorPalette(92)
prismatic::color(mycolor)
names(mycolor)=unique(data2$value)
write.table(as.data.frame(mycolor),'data/drug_color.txt',sep='\t',quote=F,row.names=F,col.names=F)
library(ggalluvial)
p <- ggplot(data2, aes(x = variable, y = link,
                       stratum = value, alluvium = flow, fill = value)) +
  geom_stratum(width = 1/8,alpha = .5) +
  geom_flow(width= 0.3,curve_type = "cubic") +
  scale_fill_manual(values = mycolor) +
  geom_text(stat = 'stratum',infer.label = TRUE,  size = 4) +
  labs(x = '', y = '') +
  theme(panel.background = element_blank(),axis.text = element_text(size=20),
        line = element_blank(), axis.text.y = element_blank(),legend.position = "none")
pdf('data/drug/Figure6A.pdf',width=18,height=22)
dev.off()

SARS3=read.csv('data/SARS_RBP_intersect.csv',header=T)
SARS3_RBP=unique(SARS3$Gene)
drug=read.csv('data/SARS_RBP_drugname.txt',sep='\t',header=T)
drug=drug[order(drug$Z_score),]

drug_target=read.csv('data/drug_humanP.txt',sep='\t',header=T)
colnames(drug_target)[2]='Target'
drug2=merge(drug,drug_target,by.x='Drug',by.y='DrugID')
drug2=drug2[order(drug2$Z_score),]

humanP_humanP=read.csv('data/humanP_humanP.txt',sep='\t',header=T)
humanP_humanP_1=humanP_humanP %>% filter(InteractionA %in% drug2$Target)
humanP_humanP_2=humanP_humanP %>% filter(InteractionB %in% drug2$Target)
humanP_humanP_2=humanP_humanP_2[,c(2,1)]
humanP_humanP2=rbind(humanP_humanP_1,humanP_humanP_2)
colnames(humanP_humanP2)=c('Target','Symbol')

sarscov=read.csv('data/SARS-CoV-2.txt',sep='\t',header=T,check.names=F)
transid=read.csv('data/SARS-CoV-2_uniprotid_transsymbol.txt',sep='\t',header=T,check.names=F)
sarscov2=merge(sarscov,transid,by.x='Gene-2',by.y='From')
sarscov2=sarscov2[,c('Gene1','To')]
colnames(sarscov2)=c('Virus','Symbol')
extrainfo=read.xlsx('data/41422_2021_581_MOESM8_ESM.xlsx',sheet=2)
extrainfo$Gene1='SARS-CoV-2'
extrainfo=extrainfo[,c(ncol(extrainfo),2)]
colnames(extrainfo)=c('Virus','Symbol')
extrainfo=extrainfo[-which(extrainfo$Symbol %in% sarscov2$Symbol),]
sarscov_all=rbind(sarscov2,extrainfo)

#Doxorubicin
drug2_test=drug2[order(drug2$Z_score),]
drug2_test2=merge(drug2_test,humanP_humanP2,by='Target')
drug2_test3=merge(drug2_test2,SARS3,by.x='Symbol',by.y='Gene')
drug2_test3=drug2_test3 %>% dplyr::select(Name,Target,Symbol,Gene1,Z_score)
colnames(drug2_test3)=c('Drug','Target','SARS_RBP','SARS_Cov','Z_score')
drug2_test3=drug2_test3[order(drug2_test3$Z_score),]
drug0='Doxorubicin'
drug2_test4=drug2_test3 %>% filter(Drug ==drug0)
drug_info=cbind(drug2_test4$Drug,'Drug')
target_info=cbind(drug2_test4$Target,'Target')
SARS_RBP_info=cbind(drug2_test4$SARS_RBP,'SARS_RBP')
SARS_Cov_info=cbind(drug2_test4$SARS_Cov,'SARS_Cov')
Nodetype=rbind(drug_info,target_info,SARS_RBP_info,SARS_Cov_info)
colnames(Nodetype)=c('Node','Type')
Nodetype=as.data.frame(Nodetype) %>% distinct()
write.table(Nodetype,paste0('data/drug/NodeType_',drug0,'.txt'),sep='\t',quote=F,row.names=F,col.names=T)
drug2_test4_1=drug2_test4[,c(1,2)] %>% distinct() %>% as.matrix()
drug2_test4_1=cbind(drug2_test4_1,'drug')
drug2_test4_2=drug2_test4[,c(2,3)] %>% distinct() %>% as.matrix()
drug2_test4_2=cbind(drug2_test4_2,'target')
drug2_test4_3=drug2_test4[,c(3,4)] %>% distinct() %>% as.matrix()
drug2_test4_3=cbind(drug2_test4_3,'sars')
drug2_test5=rbind(drug2_test4_1,drug2_test4_2,drug2_test4_3)
colnames(drug2_test5)=c('Drug','Target','SARS_Cov')
write.table(drug2_test5,paste0('data/DrugNetwork_',drug0,'.txt'),sep='\t',quote=F,row.names=F,col.names=T)

#Topotecan
drug2_test=drug2[order(drug2$Z_score),]
drug2_test2=merge(drug2_test,humanP_humanP2,by='Target')
drug2_test3=merge(drug2_test2,SARS3,by.x='Symbol',by.y='Gene')
drug2_test3=drug2_test3 %>% dplyr::select(Name,Target,Symbol,Gene1,Z_score)
colnames(drug2_test3)=c('Drug','Target','SARS_RBP','SARS_Cov','Z_score')
drug2_test3=drug2_test3[order(drug2_test3$Z_score),]
drug0='Topotecan'
drug2_test4=drug2_test3 %>% filter(Drug ==drug0)
drug_info=cbind(drug2_test4$Drug,'Drug')
target_info=cbind(drug2_test4$Target,'Target')
SARS_RBP_info=cbind(drug2_test4$SARS_RBP,'SARS_RBP')
SARS_Cov_info=cbind(drug2_test4$SARS_Cov,'SARS_Cov')
Nodetype=rbind(drug_info,target_info,SARS_RBP_info,SARS_Cov_info)
colnames(Nodetype)=c('Node','Type')
Nodetype=as.data.frame(Nodetype) %>% distinct()
write.table(Nodetype,paste0('data/NodeType_',drug0,'.txt'),sep='\t',quote=F,row.names=F,col.names=T)
drug2_test4_1=drug2_test4[,c(1,2)] %>% distinct() %>% as.matrix()
drug2_test4_1=cbind(drug2_test4_1,'Drug')
drug2_test4_2=drug2_test4[,c(2,3)] %>% distinct() %>% as.matrix()
drug2_test4_2=cbind(drug2_test4_2,'Target')
drug2_test4_3=drug2_test4[,c(3,4)] %>% distinct() %>% as.matrix()
drug2_test4_3=cbind(drug2_test4_3,'SARS_Cov')
drug2_test5=rbind(drug2_test4_1,drug2_test4_2,drug2_test4_3)
colnames(drug2_test5)=c('Drug','Target','SARS_Cov')
write.table(drug2_test5,paste0('data/DrugNetwork_',drug0,'.txt'),sep='\t',quote=F,row.names=F,col.names=T)
