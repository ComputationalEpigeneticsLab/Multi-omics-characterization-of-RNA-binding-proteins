library(openxlsx)
RBPlist=read.xlsx(paste0(datapath,'RBP.xlsx'))
RBP=RBPlist$RBP.name
RBP=RBP[!is.na(RBP)]
SARS3=read.csv('data/SARS_RBP_intersect.csv',header=T)
SARS3_RBP=unique(SARS3$Gene)

library(ggprism)
library(ggpubr)
library(ggplot2)
# BiocManager::install("ggthemes")
library(ggthemes)
library(ggprism)
library(ggpubr)
library(igraph)
library(ggsci)

humanP_humanP=read.csv('data/humanP_humanP.txt',sep='\t',header=T)
dim(humanP_humanP)
humanP_humanP_1=humanP_humanP %>% filter(InteractionA %in% RBP & InteractionB %in% RBP)
# write.table(humanP_humanP_1,'data/RBP_intersect.txt',sep='\t',quote=F,row.names=F)


SARS3=read.csv('data/SARS_RBP_intersect.csv',header=T)
SARS3_RBP=unique(SARS3$Gene)

kk=make_graph(t(humanP_humanP),directed = F)
degreemy=degree(kk,mode = "total")
closenessmy=closeness(kk,mode="in")
betweennessmy=betweenness(kk,normalized = T)

degreemy_RBP=as.data.frame(degreemy[intersect(names(degreemy),RBP)])
colnames(degreemy_RBP)='degree'
degreemy_RBP$degree=as.numeric(degreemy_RBP$degree)
degreemy_RBP$RBP=rownames(degreemy_RBP)
degreemy_RBP$SARS_interaction='Other proteins'
degreemy_RBP[which(degreemy_RBP$RBP %in% SARS3_RBP),'SARS_interaction']='SARS-CoV-2 interacting proteins'
degreemy_RBP$SARS_interaction=factor(degreemy_RBP$SARS_interaction,levels=c('SARS-CoV-2 interacting proteins','Other proteins'))
degreemy_RBP=degreemy_RBP[order(degreemy_RBP$SARS_interaction),]

closenessmy_RBP=as.data.frame(closenessmy[intersect(names(closenessmy),RBP)])
colnames(closenessmy_RBP)='closeness'
closenessmy_RBP$closeness=as.numeric(closenessmy_RBP$closeness)
closenessmy_RBP$RBP=rownames(closenessmy_RBP)
closenessmy_RBP$SARS_interaction='Other proteins'
closenessmy_RBP[which(closenessmy_RBP$RBP %in% SARS3_RBP),'SARS_interaction']='SARS-CoV-2 interacting proteins'
closenessmy_RBP$SARS_interaction=factor(closenessmy_RBP$SARS_interaction,levels=c('SARS-CoV-2 interacting proteins','Other proteins'))
closenessmy_RBP=closenessmy_RBP[order(closenessmy_RBP$SARS_interaction),]


betweennessmy_RBP=as.data.frame(betweennessmy[intersect(names(betweennessmy),RBP)])
colnames(betweennessmy_RBP)='betweenness'
betweennessmy_RBP$betweenness=as.numeric(betweennessmy_RBP$betweenness)
betweennessmy_RBP$RBP=rownames(betweennessmy_RBP)
betweennessmy_RBP$SARS_interaction='Other proteins'
betweennessmy_RBP[which(betweennessmy_RBP$RBP %in% SARS3_RBP),'SARS_interaction']='SARS-CoV-2 interacting proteins'
betweennessmy_RBP$SARS_interaction=factor(betweennessmy_RBP$SARS_interaction,levels=c('SARS-CoV-2 interacting proteins','Other proteins'))
betweennessmy_RBP=betweennessmy_RBP[order(betweennessmy_RBP$SARS_interaction),]


P1<- ggplot(degreemy_RBP, aes(x=SARS_interaction, y=log2(degree),fill=SARS_interaction)) + 
  geom_violin(trim=FALSE,color="white") + 
  geom_boxplot(width=0.2,position=position_dodge(0.9),fill='white',outlier.shape=NA)+
scale_fill_manual(values=pal_lancet('lanonc')(9)[2:1])+
ylim(5,12.5)+
theme_classic()+stat_compare_means(label = "p.format",label.x.npc = 'center')

P1

P2<- ggplot(closenessmy_RBP, aes(x=SARS_interaction, y=log2(closeness),fill=SARS_interaction)) + 
  geom_violin(trim=FALSE,color="white") +
  geom_boxplot(width=0.2,position=position_dodge(0.9),fill='white',outlier.shape=NA)+ 
 scale_fill_brewer(palette='Set1',direction=-1)+ 
scale_fill_manual(values=pal_lancet('lanonc')(9)[2:1])+

ylim(-15.6,-15.2)+
theme_classic()+stat_compare_means(label = "p.format",label.x.npc = 'center')

P2

P3<- ggplot(betweennessmy_RBP, aes(x=SARS_interaction, y=log2(betweenness),fill=SARS_interaction)) + 
  geom_violin(trim=FALSE,color="white") + 
  geom_boxplot(width=0.2,position=position_dodge(0.9),fill='white',outlier.shape=NA)+

scale_fill_manual(values=pal_lancet('lanonc')(9)[2:1])+
ylim(-20,-8)+
theme_classic()+stat_compare_means(label = "p.format",label.x.npc = 'center')

P3


library(patchwork)
pdf('Figure4_new.pdf',width=14,height=6)
P1+P2+P3+plot_layout(guides='collect',widths = c(1,1,1))
dev.off()
