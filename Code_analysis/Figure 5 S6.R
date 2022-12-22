SARS3=read.csv('data/SARS_RBP_intersect.csv',header=T)
SARS3_RBP=unique(SARS3$Gene)
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
# transform id  
map_dt <- bitr(SARS3_RBP, fromType = "SYMBOL",toType = c( "ENTREZID"),OrgDb = org.Hs.eg.db)
write.table(map_dt$ENTREZID,'data/disease299/test/RBP_SARS_id.txt',sep='\t',quote=F,row.names=F,col.names=F)

'''Linux
nohup python disease_distance.py --g1 RBP_SARS_id.txt > RBP_SARS_id.log 2>&1 &
'''

library(ggplot2)
library(tidyverse)
library(ragg)
data1<-read.csv("data/disease299/RBP_SARS_id_separation_results.txt",stringsAsFactors = F,header = F,sep="\t")
colnames(data1)=c('mygeneset','mygeneset_len','disease_name','disease_geneset_len',
                  'd_A','d_B','d_AB','s_AB')
data1$mygeneset=gsub('.txt','',data1$mygeneset)
data1$disease_name=gsub('.txt','',data1$disease_name)
data1=data1 %>% dplyr::select('disease_name','disease_geneset_len','s_AB')
type<-read.csv("data/disease299/disease_type.txt",stringsAsFactors = F,header = T,sep="\t")
data<-merge(data1,type,by="disease_name")
colnames(data)<-c("individual","number","value","group")
datagroup <- data$group %>% unique()
allplotdata <- tibble('group' = datagroup,
                      'individual' = paste0('empty_individual_', seq_along(datagroup)),
                      'value' = 0,'number'=0) %>% 
  bind_rows(data) %>% arrange(group) %>% mutate(xid = 1:n()) %>% 
  mutate(angle = 90 - 360 * (xid - 0.5) / n()) %>% 
  mutate(hjust = ifelse(angle < -90, 1, 0)) %>% 
  mutate(angle = ifelse(angle < -90, angle+180, angle)) 

disease_group_col=c("#D06862","#8CC4E5","#778CB7","#D17EAE",
                    "#A2D6B9","#637B54","#FED13F","#9E66A5",
                    "#F79E6D","#63CCE3")
names(disease_group_col)=sort(unique(re$group))
disease_group_col

firstxid <- which(str_detect(allplotdata$individual, pattern = "empty_individual")) # 1 12 43 58
re<-allplotdata[-firstxid,]
segment_data <- data.frame('from' = firstxid + 1,
                           'to' = c(c(firstxid - 1)[-1], nrow(allplotdata)),
                           'label' = datagroup) %>% 
  mutate(labelx = as.integer((from + to)/2))
re$value=as.numeric(re$value)
allplotdata2=allplotdata %>% filter(value < 0)

p=ggplot(re,aes(x = xid, y = value),colour=group) + 
  geom_point(data = allplotdata, 
             aes(x = xid, y = value,size=number,color=group),
             stat = 'identity',alpha=1,shape=20)+
  geom_point(data = allplotdata2, 
             aes(x = xid, y = value,size=number),
             stat = 'identity',alpha=1,shape=21)+
  theme_minimal()+
  scale_x_continuous(expand = c(0, 12)) +
  scale_y_continuous(breaks = c(-0.6,-0.4,-0.2,0,0.1,0.2,0.4,0.6,0.8),limits = c(-0.6,0.8))+
  scale_color_manual(values = c("#D06862","#8CC4E5","#778CB7","#D17EAE",
                                "#A2D6B9","#637B54","#FED13F","#9E66A5",
                                "#F79E6D","#63CCE3"))+
  scale_size(range = c(2,10))
p=p+coord_polar(theta = "x", start =-pi/1.75,direction = -1)

library(ggrepel)
mygeneall=as.data.frame(re[which(re$value < 0 ),'individual'])
mygeneall=as.character(mygeneall[,1])
mygeneall
p=p+geom_text_repel(
  data = re[which(re$individual %in% mygeneall),],
  aes(label = individual,color=group),
  size =3,
  segment.color = "black", show.legend = FALSE )+
  coord_polar(theta = "x", start =-pi/1.75,direction = -1)

print(disease_group_col)
pdf('disease299_point.pdf',width=8,height=8)
print(p)
dev.off()


#kidney
data1<-read.csv("RBP_SARS_id_separation_results.txt",stringsAsFactors = F,header = F,sep="\t")
colnames(data1)=c('mygeneset','mygeneset_len','disease_name','disease_geneset_len',
                  'd_A','d_B','d_AB','s_AB')
data1=data1[order(data1$s_AB),]
disease0=data1$disease_name[1]
disease_geneset=read.csv(paste0('data/disease299/disease299/',disease0),sep='\t',header=F)
disease_geneset2=as.numeric(disease_geneset[,1])
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(stringr)
map_dt <- bitr(disease_geneset2, fromType = "ENTREZID",toType = c( "SYMBOL"),OrgDb = org.Hs.eg.db)
map_dt$Disease=gsub('.txt','',str_to_title(disease0))
colnames(map_dt)=c('Disease_ENTREZID','Disease_SYMBOL','Disease')

SARS3=read.csv('data/SARS_RBP_intersect.csv',header=T)
SARS3_RBP=unique(SARS3$Gene)
humanP_humanP=read.csv('data/humanP_humanP.txt',sep='\t',header=T)
humanP_humanP_1=humanP_humanP %>% filter(InteractionA %in% map_dt$Disease_SYMBOL)
humanP_humanP_2=humanP_humanP %>% filter(InteractionB %in% map_dt$Disease_SYMBOL)
humanP_humanP_2=humanP_humanP_2[,c(2,1)]
humanP_humanP2=rbind(humanP_humanP_1,humanP_humanP_2)
colnames(humanP_humanP2)=c('Disease_SYMBOL','SYMBOL')
humanP_humanP2=as.data.frame(humanP_humanP2) %>% distinct()
humanP_humanP2$Class=apply(humanP_humanP2,1,function(x){
  a=sort(c(x[1],x[2]))
  b=paste(a,collapse='_')
  return(b)
})
humanP_humanP3=humanP_humanP2[!duplicated(humanP_humanP2$Class),]
humanP_humanP3=humanP_humanP3[,1:2]
humanP_humanP3=humanP_humanP3 %>% filter(SYMBOL %in% SARS3_RBP)
map_dt2=merge(map_dt,humanP_humanP3,by='Disease_SYMBOL')

SARS_RBP_intersect=read.csv('data/SARS_RBP_intersect.csv',header=T)
SARS_RBP_intersect=SARS_RBP_intersect %>% dplyr::select(Gene,Gene1) %>% distinct()
colnames(SARS_RBP_intersect)=c('RBP','SARS_Cov')

map_dt3=merge(map_dt2,SARS_RBP_intersect,by.x='SYMBOL',by.y='RBP',all.x=T)
drug_info=cbind(map_dt3$Disease_SYMBOL,'Disease_SYMBOL')
target_info=cbind(map_dt3$Disease,'Disease')
SARS_RBP_info=cbind(map_dt3$SYMBOL,'SARS_RBP')
SARS_info=cbind(map_dt3$SARS_Cov,'SARS-Cov-2')
Nodetype=rbind(drug_info,target_info,SARS_RBP_info,SARS_info)
colnames(Nodetype)=c('Node','Type')
Nodetype=as.data.frame(Nodetype) %>% distinct()
Nodetype[which(Nodetype$Node %in% map_dt$Disease_SYMBOL),'Type']='Disease_SYMBOL'
Nodetype=Nodetype %>% distinct()
Nodetype[which(Nodetype$Node %in% intersect(map_dt$Disease_SYMBOL,SARS3_RBP)),'Type']='Double_SYMBOL'
write.table(Nodetype,paste0('data/disease299/NodeType_',gsub('.txt','',str_to_title(disease0)),'_1.txt'),sep='\t',quote=F,row.names=F,col.names=T)

map_dt3_1=map_dt3[,c('Disease_SYMBOL','Disease')] %>% distinct() %>% as.matrix()
map_dt3_1=cbind(map_dt3_1,'Disease')

map_dt3_2=map_dt3[,c('Disease_SYMBOL','SYMBOL')] %>% distinct() %>% as.matrix()
map_dt3_2=cbind(map_dt3_2,'SYMBOL')

map_dt3_3=SARS_RBP_intersect[which(SARS_RBP_intersect$RBP %in% c(map_dt3$SYMBOL,map_dt3$Disease_SYMBOL)),] %>% distinct() %>% as.matrix()
map_dt3_3=cbind(map_dt3_3,'SARS-Cov-2')

map_dt4=rbind(map_dt3_1,map_dt3_2,map_dt3_3)
colnames(map_dt4)=c('Node1','Node2','Edge')
write.table(map_dt4,paste0('Network_',gsub('.txt','',str_to_title(disease0)),'_1.txt'),sep='\t',quote=F,row.names=F,col.names=T)
#these file were submitted to Cytoscape for visualization

#esophagus
data1<-read.csv("RBP_SARS_id_separation_results.txt",stringsAsFactors = F,header = F,sep="\t")
colnames(data1)=c('mygeneset','mygeneset_len','disease_name','disease_geneset_len',
                  'd_A','d_B','d_AB','s_AB')
data1=data1[order(data1$s_AB),]
disease0=data1$disease_name[6]
disease_geneset=read.csv(paste0('disease299/',disease0),sep='\t',header=F)
disease_geneset2=as.numeric(disease_geneset[,1])
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(stringr)
# transform id  
map_dt <- bitr(disease_geneset2, fromType = "ENTREZID",toType = c( "SYMBOL"),OrgDb = org.Hs.eg.db)
map_dt$Disease=gsub('.txt','',str_to_title(disease0))
colnames(map_dt)=c('Disease_ENTREZID','Disease_SYMBOL','Disease')

SARS3=read.csv('data/SARS_RBP_intersect.csv',header=T)
SARS3_RBP=unique(SARS3$Gene)
humanP_humanP=read.csv('data/humanP_humanP.txt',sep='\t',header=T)
humanP_humanP_1=humanP_humanP %>% filter(InteractionA %in% map_dt$Disease_SYMBOL)
humanP_humanP_2=humanP_humanP %>% filter(InteractionB %in% map_dt$Disease_SYMBOL)
humanP_humanP_2=humanP_humanP_2[,c(2,1)]
humanP_humanP2=rbind(humanP_humanP_1,humanP_humanP_2)
colnames(humanP_humanP2)=c('Disease_SYMBOL','SYMBOL')
humanP_humanP2=as.data.frame(humanP_humanP2) %>% distinct()
humanP_humanP2$Class=apply(humanP_humanP2,1,function(x){
  a=sort(c(x[1],x[2]))
  b=paste(a,collapse='_')
  return(b)
})
humanP_humanP3=humanP_humanP2[!duplicated(humanP_humanP2$Class),]
humanP_humanP3=humanP_humanP3[,1:2]
humanP_humanP3=humanP_humanP3 %>% filter(SYMBOL %in% SARS3_RBP)
map_dt2=merge(map_dt,humanP_humanP3,by='Disease_SYMBOL')

SARS_RBP_intersect=read.csv('data/SARS_RBP_intersect.csv',header=T)
SARS_RBP_intersect=SARS_RBP_intersect %>% dplyr::select(Gene,Gene1) %>% distinct()
colnames(SARS_RBP_intersect)=c('RBP','SARS_Cov')

map_dt3=merge(map_dt2,SARS_RBP_intersect,by.x='SYMBOL',by.y='RBP',all.x=T)
drug_info=cbind(map_dt3$Disease_SYMBOL,'Disease_SYMBOL')
target_info=cbind(map_dt3$Disease,'Disease')
SARS_RBP_info=cbind(map_dt3$SYMBOL,'SARS_RBP')
SARS_info=cbind(map_dt3$SARS_Cov,'SARS-Cov-2')
Nodetype=rbind(drug_info,target_info,SARS_RBP_info,SARS_info)
colnames(Nodetype)=c('Node','Type')
Nodetype=as.data.frame(Nodetype) %>% distinct()
Nodetype[which(Nodetype$Node %in% map_dt$Disease_SYMBOL),'Type']='Disease_SYMBOL'
Nodetype=Nodetype %>% distinct()
Nodetype[which(Nodetype$Node %in% intersect(map_dt$Disease_SYMBOL,SARS3_RBP)),'Type']='Double_SYMBOL'
write.table(Nodetype,paste0('NodeType_',gsub('.txt','',str_to_title(disease0)),'_1.txt'),sep='\t',quote=F,row.names=F,col.names=T)

map_dt3_1=map_dt3[,c('Disease_SYMBOL','Disease')] %>% distinct() %>% as.matrix()
map_dt3_1=cbind(map_dt3_1,'Disease')

map_dt3_2=map_dt3[,c('Disease_SYMBOL','SYMBOL')] %>% distinct() %>% as.matrix()
map_dt3_2=cbind(map_dt3_2,'SYMBOL')

map_dt3_3=SARS_RBP_intersect[which(SARS_RBP_intersect$RBP %in% c(map_dt3$SYMBOL,map_dt3$Disease_SYMBOL)),] %>% distinct() %>% as.matrix()
map_dt3_3=cbind(map_dt3_3,'SARS-Cov-2')

map_dt4=rbind(map_dt3_1,map_dt3_2,map_dt3_3)
colnames(map_dt4)=c('Node1','Node2','Edge')
write.table(map_dt4,paste0('Network_',gsub('.txt','',str_to_title(disease0)),'_1.txt'),sep='\t',quote=F,row.names=F,col.names=T)
#these file were submitted to Cytoscape for visualization

#urinary1
data1<-read.csv("RBP_SARS_id_separation_results.txt",stringsAsFactors = F,header = F,sep="\t")
colnames(data1)=c('mygeneset','mygeneset_len','disease_name','disease_geneset_len',
                  'd_A','d_B','d_AB','s_AB')
data1=data1[order(data1$s_AB),]
disease0=data1$disease_name[5]
disease_geneset=read.csv(paste0('disease299/',disease0),sep='\t',header=F)
disease_geneset2=as.numeric(disease_geneset[,1])
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(stringr)
map_dt <- bitr(disease_geneset2, fromType = "ENTREZID",toType = c( "SYMBOL"),OrgDb = org.Hs.eg.db)
map_dt$Disease=gsub('.txt','',str_to_title(disease0))
colnames(map_dt)=c('Disease_ENTREZID','Disease_SYMBOL','Disease')

SARS3=read.csv('data/SARS_RBP_intersect.csv',header=T)
SARS3_RBP=unique(SARS3$Gene)
humanP_humanP=read.csv('humanP_humanP.txt',sep='\t',header=T)
humanP_humanP_1=humanP_humanP %>% filter(InteractionA %in% map_dt$Disease_SYMBOL)
humanP_humanP_2=humanP_humanP %>% filter(InteractionB %in% map_dt$Disease_SYMBOL)
humanP_humanP_2=humanP_humanP_2[,c(2,1)]
humanP_humanP2=rbind(humanP_humanP_1,humanP_humanP_2)
colnames(humanP_humanP2)=c('Disease_SYMBOL','SYMBOL')
humanP_humanP2=as.data.frame(humanP_humanP2) %>% distinct()
humanP_humanP2$Class=apply(humanP_humanP2,1,function(x){
  a=sort(c(x[1],x[2]))
  b=paste(a,collapse='_')
  return(b)
})
humanP_humanP3=humanP_humanP2[!duplicated(humanP_humanP2$Class),]
humanP_humanP3=humanP_humanP3[,1:2]
humanP_humanP3=humanP_humanP3 %>% filter(SYMBOL %in% SARS3_RBP)
map_dt2=merge(map_dt,humanP_humanP3,by='Disease_SYMBOL')

SARS_RBP_intersect=read.csv('SARS_RBP_intersect.csv',header=T)
SARS_RBP_intersect=SARS_RBP_intersect %>% dplyr::select(Gene,Gene1) %>% distinct()
colnames(SARS_RBP_intersect)=c('RBP','SARS_Cov')

map_dt3=merge(map_dt2,SARS_RBP_intersect,by.x='SYMBOL',by.y='RBP',all.x=T)
drug_info=cbind(map_dt3$Disease_SYMBOL,'Disease_SYMBOL')
target_info=cbind(map_dt3$Disease,'Disease')
SARS_RBP_info=cbind(map_dt3$SYMBOL,'SARS_RBP')
SARS_info=cbind(map_dt3$SARS_Cov,'SARS-Cov-2')
Nodetype=rbind(drug_info,target_info,SARS_RBP_info,SARS_info)
colnames(Nodetype)=c('Node','Type')
Nodetype=as.data.frame(Nodetype) %>% distinct()
Nodetype[which(Nodetype$Node %in% map_dt$Disease_SYMBOL),'Type']='Disease_SYMBOL'
Nodetype=Nodetype %>% distinct()
Nodetype[which(Nodetype$Node %in% intersect(map_dt$Disease_SYMBOL,SARS3_RBP)),'Type']='Double_SYMBOL'
write.table(Nodetype,paste0('NodeType_',gsub('.txt','',str_to_title(disease0)),'_1.txt'),sep='\t',quote=F,row.names=F,col.names=T)

map_dt3_1=map_dt3[,c('Disease_SYMBOL','Disease')] %>% distinct() %>% as.matrix()
map_dt3_1=cbind(map_dt3_1,'Disease')

map_dt3_2=map_dt3[,c('Disease_SYMBOL','SYMBOL')] %>% distinct() %>% as.matrix()
map_dt3_2=cbind(map_dt3_2,'SYMBOL')

map_dt3_3=SARS_RBP_intersect[which(SARS_RBP_intersect$RBP %in% c(map_dt3$SYMBOL,map_dt3$Disease_SYMBOL)),] %>% distinct() %>% as.matrix()
map_dt3_3=cbind(map_dt3_3,'SARS-Cov-2')

map_dt4=rbind(map_dt3_1,map_dt3_2,map_dt3_3)
colnames(map_dt4)=c('Node1','Node2','Edge')
write.table(map_dt4,paste0('data/disease299/Network_',gsub('.txt','',str_to_title(disease0)),'_1.txt'),sep='\t',quote=F,row.names=F,col.names=T)
map_dt3_1=map_dt3[,c('Disease_SYMBOL','Disease')] %>% distinct() %>% as.matrix()
map_dt3_1=cbind(map_dt3_1,'Disease')

map_dt3_2=map_dt3[,c('Disease_SYMBOL','SYMBOL')] %>% distinct() %>% as.matrix()
map_dt3_2=cbind(map_dt3_2,'SYMBOL')

map_dt3_3=SARS_RBP_intersect[which(SARS_RBP_intersect$RBP %in% c(map_dt3$SYMBOL,map_dt3$Disease_SYMBOL)),] %>% distinct() %>% as.matrix()
map_dt3_3=cbind(map_dt3_3,'SARS-Cov-2')

map_dt4=rbind(map_dt3_1,map_dt3_2,map_dt3_3)
colnames(map_dt4)=c('Node1','Node2','Edge')

#urinary2
data1<-read.csv("RBP_SARS_id_separation_results.txt",stringsAsFactors = F,header = F,sep="\t")
colnames(data1)=c('mygeneset','mygeneset_len','disease_name','disease_geneset_len',
                  'd_A','d_B','d_AB','s_AB')
data1=data1[order(data1$s_AB),]
disease0=data1$disease_name[7]
disease_geneset=read.csv(paste0('data/',disease0),sep='\t',header=F)
disease_geneset2=as.numeric(disease_geneset[,1])
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(stringr)
map_dt <- bitr(disease_geneset2, fromType = "ENTREZID",toType = c( "SYMBOL"),OrgDb = org.Hs.eg.db)
map_dt$Disease=gsub('.txt','',str_to_title(disease0))
colnames(map_dt)=c('Disease_ENTREZID','Disease_SYMBOL','Disease')

SARS3=read.csv('data/SARS_RBP_intersect.csv',header=T)
SARS3_RBP=unique(SARS3$Gene)
humanP_humanP=read.csv('data/humanP_humanP.txt',sep='\t',header=T)
humanP_humanP_1=humanP_humanP %>% filter(InteractionA %in% map_dt$Disease_SYMBOL)
humanP_humanP_2=humanP_humanP %>% filter(InteractionB %in% map_dt$Disease_SYMBOL)
humanP_humanP_2=humanP_humanP_2[,c(2,1)]
humanP_humanP2=rbind(humanP_humanP_1,humanP_humanP_2)
colnames(humanP_humanP2)=c('Disease_SYMBOL','SYMBOL')
humanP_humanP2=as.data.frame(humanP_humanP2) %>% distinct()
humanP_humanP2$Class=apply(humanP_humanP2,1,function(x){
  a=sort(c(x[1],x[2]))
  b=paste(a,collapse='_')
  return(b)
})
humanP_humanP3=humanP_humanP2[!duplicated(humanP_humanP2$Class),]
humanP_humanP3=humanP_humanP3[,1:2]
humanP_humanP3=humanP_humanP3 %>% filter(SYMBOL %in% SARS3_RBP)
map_dt2=merge(map_dt,humanP_humanP3,by='Disease_SYMBOL')

SARS_RBP_intersect=read.csv('data/SARS_RBP_intersect.csv',header=T)
SARS_RBP_intersect=SARS_RBP_intersect %>% dplyr::select(Gene,Gene1) %>% distinct()
colnames(SARS_RBP_intersect)=c('RBP','SARS_Cov')

map_dt3=merge(map_dt2,SARS_RBP_intersect,by.x='SYMBOL',by.y='RBP',all.x=T)
drug_info=cbind(map_dt3$Disease_SYMBOL,'Disease_SYMBOL')
target_info=cbind(map_dt3$Disease,'Disease')
SARS_RBP_info=cbind(map_dt3$SYMBOL,'SARS_RBP')
SARS_info=cbind(map_dt3$SARS_Cov,'SARS-Cov-2')
Nodetype=rbind(drug_info,target_info,SARS_RBP_info,SARS_info)
colnames(Nodetype)=c('Node','Type')
Nodetype=as.data.frame(Nodetype) %>% distinct()
Nodetype[which(Nodetype$Node %in% map_dt$Disease_SYMBOL),'Type']='Disease_SYMBOL'
Nodetype=Nodetype %>% distinct()
Nodetype[which(Nodetype$Node %in% intersect(map_dt$Disease_SYMBOL,SARS3_RBP)),'Type']='Double_SYMBOL'
write.table(Nodetype,paste0('data/NodeType_',gsub('.txt','',str_to_title(disease0)),'_1.txt'),sep='\t',quote=F,row.names=F,col.names=T)

map_dt3_1=map_dt3[,c('Disease_SYMBOL','Disease')] %>% distinct() %>% as.matrix()
map_dt3_1=cbind(map_dt3_1,'Disease')

map_dt3_2=map_dt3[,c('Disease_SYMBOL','SYMBOL')] %>% distinct() %>% as.matrix()
map_dt3_2=cbind(map_dt3_2,'SYMBOL')

map_dt3_3=SARS_RBP_intersect[which(SARS_RBP_intersect$RBP %in% c(map_dt3$SYMBOL,map_dt3$Disease_SYMBOL)),] %>% distinct() %>% as.matrix()
map_dt3_3=cbind(map_dt3_3,'SARS-Cov-2')

map_dt4=rbind(map_dt3_1,map_dt3_2,map_dt3_3)
colnames(map_dt4)=c('Node1','Node2','Edge')
write.table(map_dt4,paste0('data/Network_',gsub('.txt','',str_to_title(disease0)),'_1.txt'),sep='\t',quote=F,row.names=F,col.names=T)
map_dt3_1=map_dt3[,c('Disease_SYMBOL','Disease')] %>% distinct() %>% as.matrix()
map_dt3_1=cbind(map_dt3_1,'Disease')

map_dt3_2=map_dt3[,c('Disease_SYMBOL','SYMBOL')] %>% distinct() %>% as.matrix()
map_dt3_2=cbind(map_dt3_2,'SYMBOL')

map_dt3_3=SARS_RBP_intersect[which(SARS_RBP_intersect$RBP %in% c(map_dt3$SYMBOL,map_dt3$Disease_SYMBOL)),] %>% distinct() %>% as.matrix()
map_dt3_3=cbind(map_dt3_3,'SARS-Cov-2')

map_dt4=rbind(map_dt3_1,map_dt3_2,map_dt3_3)
colnames(map_dt4)=c('Node1','Node2','Edge')
write.table(map_dt4,paste0('data/Network_',gsub('.txt','',str_to_title(disease0)),'_1.txt'),sep='\t',quote=F,row.names=F,col.names=T)
write.table(map_dt4,paste0('Network_',gsub('.txt','',str_to_title(disease0)),'_1.txt'),sep='\t',quote=F,row.names=F,col.names=T)
#these file were submitted to Cytoscape for visualization

#leukaemia
data1<-read.csv("data/RBP_SARS_id_separation_results.txt",stringsAsFactors = F,header = F,sep="\t")
colnames(data1)=c('mygeneset','mygeneset_len','disease_name','disease_geneset_len',
                  'd_A','d_B','d_AB','s_AB')
data1=data1[order(data1$s_AB),]
disease0=data1$disease_name[10]
disease_geneset=read.csv(paste0('data/',disease0),sep='\t',header=F)
disease_geneset2=as.numeric(disease_geneset[,1])
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(stringr)
# transform id  
map_dt <- bitr(disease_geneset2, fromType = "ENTREZID",toType = c( "SYMBOL"),OrgDb = org.Hs.eg.db)
map_dt$Disease=gsub('.txt','',str_to_title(disease0))
colnames(map_dt)=c('Disease_ENTREZID','Disease_SYMBOL','Disease')

SARS3=read.csv('data/SARS_RBP_intersect.csv',header=T)
SARS3_RBP=unique(SARS3$Gene)
humanP_humanP=read.csv('data/humanP_humanP.txt',sep='\t',header=T)
humanP_humanP_1=humanP_humanP %>% filter(InteractionA %in% map_dt$Disease_SYMBOL)
humanP_humanP_2=humanP_humanP %>% filter(InteractionB %in% map_dt$Disease_SYMBOL)
humanP_humanP_2=humanP_humanP_2[,c(2,1)]
humanP_humanP2=rbind(humanP_humanP_1,humanP_humanP_2)
colnames(humanP_humanP2)=c('Disease_SYMBOL','SYMBOL')
humanP_humanP2=as.data.frame(humanP_humanP2) %>% distinct()
humanP_humanP2$Class=apply(humanP_humanP2,1,function(x){
  a=sort(c(x[1],x[2]))
  b=paste(a,collapse='_')
  return(b)
})
humanP_humanP3=humanP_humanP2[!duplicated(humanP_humanP2$Class),]
humanP_humanP3=humanP_humanP3[,1:2]
humanP_humanP3=humanP_humanP3 %>% filter(SYMBOL %in% SARS3_RBP)
map_dt2=merge(map_dt,humanP_humanP3,by='Disease_SYMBOL')

SARS_RBP_intersect=read.csv('data/SARS_RBP_intersect.csv',header=T)
SARS_RBP_intersect=SARS_RBP_intersect %>% dplyr::select(Gene,Gene1) %>% distinct()
colnames(SARS_RBP_intersect)=c('RBP','SARS_Cov')

map_dt3=merge(map_dt2,SARS_RBP_intersect,by.x='SYMBOL',by.y='RBP',all.x=T)
drug_info=cbind(map_dt3$Disease_SYMBOL,'Disease_SYMBOL')
target_info=cbind(map_dt3$Disease,'Disease')
SARS_RBP_info=cbind(map_dt3$SYMBOL,'SARS_RBP')
SARS_info=cbind(map_dt3$SARS_Cov,'SARS-Cov-2')
Nodetype=rbind(drug_info,target_info,SARS_RBP_info,SARS_info)
colnames(Nodetype)=c('Node','Type')
Nodetype=as.data.frame(Nodetype) %>% distinct()
Nodetype[which(Nodetype$Node %in% map_dt$Disease_SYMBOL),'Type']='Disease_SYMBOL'
Nodetype=Nodetype %>% distinct()
Nodetype[which(Nodetype$Node %in% intersect(map_dt$Disease_SYMBOL,SARS3_RBP)),'Type']='Double_SYMBOL'
write.table(Nodetype,paste0('data/NodeType_',gsub('.txt','',str_to_title(disease0)),'_1.txt'),sep='\t',quote=F,row.names=F,col.names=T)

map_dt3_1=map_dt3[,c('Disease_SYMBOL','Disease')] %>% distinct() %>% as.matrix()
map_dt3_1=cbind(map_dt3_1,'Disease')

map_dt3_2=map_dt3[,c('Disease_SYMBOL','SYMBOL')] %>% distinct() %>% as.matrix()
map_dt3_2=cbind(map_dt3_2,'SYMBOL')

map_dt3_3=SARS_RBP_intersect[which(SARS_RBP_intersect$RBP %in% c(map_dt3$SYMBOL,map_dt3$Disease_SYMBOL)),] %>% distinct() %>% as.matrix()
map_dt3_3=cbind(map_dt3_3,'SARS-Cov-2')

map_dt4=rbind(map_dt3_1,map_dt3_2,map_dt3_3)
colnames(map_dt4)=c('Node1','Node2','Edge')
write.table(map_dt4,paste0('data/Network_',gsub('.txt','',str_to_title(disease0)),'_1.txt'),sep='\t',quote=F,row.names=F,col.names=T)
map_dt3_1=map_dt3[,c('Disease_SYMBOL','Disease')] %>% distinct() %>% as.matrix()
map_dt3_1=cbind(map_dt3_1,'Disease')

map_dt3_2=map_dt3[,c('Disease_SYMBOL','SYMBOL')] %>% distinct() %>% as.matrix()
map_dt3_2=cbind(map_dt3_2,'SYMBOL')

map_dt3_3=SARS_RBP_intersect[which(SARS_RBP_intersect$RBP %in% c(map_dt3$SYMBOL,map_dt3$Disease_SYMBOL)),] %>% distinct() %>% as.matrix()
map_dt3_3=cbind(map_dt3_3,'SARS-Cov-2')

map_dt4=rbind(map_dt3_1,map_dt3_2,map_dt3_3)
colnames(map_dt4)=c('Node1','Node2','Edge')
write.table(map_dt4,paste0('data/Network_',gsub('.txt','',str_to_title(disease0)),'_1.txt'),sep='\t',quote=F,row.names=F,col.names=T)
#these file were submitted to Cytoscape for visualization