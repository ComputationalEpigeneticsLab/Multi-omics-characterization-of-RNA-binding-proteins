expdata2=reshape2::melt(as.matrix(expdata))
colnames(expdata2)=c('Gene','PatientID','Expression')
expdata2=merge(expdata2,sampletype,by='PatientID')
Disease_development= pal_startrek('uniform')(4)[c(3,2,4,1)]

Geneselect=c('LARP1','PA2G4')
pall=list()
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
  p=ggplot(data=df2,aes(x=Status,y=log2(Expression),color=Status))+#Sample_Type
    stat_boxplot(geom="errorbar",width=0.15,position = position_dodge(0.9))+
    geom_boxplot(position = position_dodge(0.9))+
    geom_jitter(aes(color=Status),size=1)+
    labs(title = paste0(Geneselect0,'_mRNA'))+xlab(NULL)+
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
  pall[[i]] = p}

library(patchwork)
pdf('FigureS3.pdf')
pall[[1]]+pall[[2]]+plot_layout(ncol=2,guides='collect')
dev.off()