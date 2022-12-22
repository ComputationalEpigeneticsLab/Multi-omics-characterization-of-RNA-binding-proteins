library(ggprism)
library(ggpubr)
library(ggplot2)
library(ggthemes)
library(ggprism)
library(ggpubr)
library(igraph)
library(ggsci)

my_theme2<-theme_prism()+
  theme(axis.text.x=element_text(colour="black",size=10),
        axis.title.x=element_text(size = 13),
        axis.text.y=element_text(size=10,face="plain"),
        axis.title.y=element_text(size = 13,),
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=1),
        legend.text=element_text(  colour="black",
                                   size=12),
        legend.title=element_text( colour="black",
                                   size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

SARS3=read.csv('data/SARS_RBP_intersect.csv',header=T)
SARS3_RBP=unique(SARS3$Gene)


#HumanNet-FN
humanP_humanP=read.csv('data/HumanNet-FN_symbol.tsv',sep='\t',header=T)
dim(humanP_humanP)
humanP_humanP=humanP_humanP[,4:5]

kk=make_graph(t(humanP_humanP),directed = F)#graph_from_data_frame(humanP_humanP)
degreemy=degree(kk,mode = "total")
# degreemy_distribution=degree.distribution(kk) 
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
  scale_fill_manual(values = c("#3A5FCD","#EE2C2C"))+
  theme_classic()+stat_compare_means()

P2<- ggplot(closenessmy_RBP, aes(x=SARS_interaction, y=log2(closeness),fill=SARS_interaction)) + 
  geom_violin(trim=FALSE,color="white") +
  geom_boxplot(width=0.2,position=position_dodge(0.9),fill='white',outlier.shape=NA)+``
  scale_fill_manual(values = c("#3A5FCD","#EE2C2C"))+ 
  theme_classic()+stat_compare_means()

P3<- ggplot(betweennessmy_RBP, aes(x=SARS_interaction, y=log2(betweenness),fill=SARS_interaction)) + 
  geom_violin(trim=FALSE,color="white") +
  geom_boxplot(width=0.2,position=position_dodge(0.9),fill='white',outlier.shape=NA)+ 
  scale_fill_manual(values = c("#3A5FCD","#EE2C2C"))+
  theme_classic()+stat_compare_means()


#HumanNet-XC
humanP_humanP=read.csv('data/HumanNet-XC_symbol.tsv',sep='\t',header=T)
dim(humanP_humanP)
humanP_humanP=humanP_humanP[,4:5]

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


P4<- ggplot(degreemy_RBP, aes(x=SARS_interaction, y=log2(degree),fill=SARS_interaction)) + 
  geom_violin(trim=FALSE,color="white") +
  geom_boxplot(width=0.2,position=position_dodge(0.9),fill='white',outlier.shape=NA)+
  scale_fill_manual(values = c("#3A5FCD","#EE2C2C"))+
  theme_classic()+stat_compare_means()

P5<- ggplot(closenessmy_RBP, aes(x=SARS_interaction, y=log2(closeness),fill=SARS_interaction)) + 
  geom_violin(trim=FALSE,color="white") +
  geom_boxplot(width=0.2,position=position_dodge(0.9),fill='white',outlier.shape=NA)+
  scale_fill_manual(values = c("#3A5FCD","#EE2C2C"))+
  theme_classic()+stat_compare_means()

P6<- ggplot(betweennessmy_RBP, aes(x=SARS_interaction, y=log2(betweenness),fill=SARS_interaction)) + 
  geom_violin(trim=FALSE,color="white") +
  geom_boxplot(width=0.2,position=position_dodge(0.9),fill='white',outlier.shape=NA)+
  scale_fill_manual(values = c("#3A5FCD","#EE2C2C"))+
  theme_classic()+stat_compare_means()

save.image('XC_symbol_igraph.RData')


#HS-PI
humanP_humanP=read.csv('data/HS-PI_symbol.tsv',sep='\t',header=T)
dim(humanP_humanP)
humanP_humanP=humanP_humanP[,4:5]

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



P7<- ggplot(degreemy_RBP, aes(x=SARS_interaction, y=log2(degree),fill=SARS_interaction)) + 
  geom_violin(trim=FALSE,color="white") +
  geom_boxplot(width=0.2,position=position_dodge(0.9),fill='white',outlier.shape=NA)+
  scale_fill_manual(values = c("#3A5FCD","#EE2C2C"))+
  theme_classic()+stat_compare_means()

P8<- ggplot(closenessmy_RBP, aes(x=SARS_interaction, y=log2(closeness),fill=SARS_interaction)) + 
  geom_violin(trim=FALSE,color="white") +
  geom_boxplot(width=0.2,position=position_dodge(0.9),fill='white',outlier.shape=NA)+
  scale_fill_manual(values = c("#3A5FCD","#EE2C2C"))+
  theme_classic()+stat_compare_means()

P9<- ggplot(betweennessmy_RBP, aes(x=SARS_interaction, y=log2(betweenness),fill=SARS_interaction)) + 
  geom_violin(trim=FALSE,color="white") +
  geom_boxplot(width=0.2,position=position_dodge(0.9),fill='white',outlier.shape=NA)+
  scale_fill_manual(values = c("#3A5FCD","#EE2C2C"))+
  theme_classic()+stat_compare_means()

library(patchwork)
pdf('S4.pdf',width=10,height=8)
P1+P2+P3+P4+P5+P6+P7+P8+P9+plot_layout(ncol=3) & theme(legend.position='none')
dev.off()