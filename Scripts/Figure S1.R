compare=c('Mild-Asymptomatic',
          'Severe-Asymptomatic',
          'Critical-Asymptomatic',
          'Severe-Mild',
          'Critical-Mild',
          'Critical-Severe')


#mRNA
plotall=list()

for (i in 1:length(compare)) {
  compare0=compare[i]
  res2=read.csv(paste0('data/DEresult/mRNA_',compare0,'_DE.txt'),sep='\t',header=T)
  res2$adj.P.Val=as.numeric(res2$adj.P.Val)
  res2$logFC=as.numeric(res2$logFC)
  threshold_new_res2 <- as.factor(ifelse(res2$adj.P.Val < 0.05 &
                                           abs(res2$logFC) > 0  ,
                                         ifelse(res2$logFC > 0 ,'Up','Down'),'Not'))
  res2$threshold=as.character(threshold_new_res2)
  alldown=length(which(res2$threshold == 'Down'))
  allup=length(which(res2$threshold == 'Up'))
  res2$Size=1
  res2[which(rownames(res2) %in% RBP & res2$threshold != 'Not'),'Size'] =2
  res2$alpha=0.5
  res2[which(rownames(res2) %in% RBP & res2$threshold != 'Not'),'alpha'] =0.8
  res2$shape='Other'
  res2[which(rownames(res2) %in% RBP & res2$threshold != 'Not'),'shape'] ='RBP'
  RBPdown=length(which(res2$threshold == 'Down' & res2$shape == 'RBP'))
  RBPup=length(which(res2$threshold == 'Up' & res2$shape == 'RBP'))
  
  ymax=max(-log10(res2$adj.P.Val))
  
  p=ggplot(res2,aes(x=logFC,y=-log10(adj.P.Val),colour=threshold)) +
    xlab("log2(Fold Change)")+ylab("-log10(padjust)") +ggtitle(compare0)+
    geom_point(aes(alpha=alpha,size = Size,shape = shape))+
    scale_size_continuous(range = c(2,3))+
    xlim(-4,4) +
    scale_color_manual(values =c("#0072B5","grey","#BC3C28"))+
    geom_vline(xintercept = c(0), lty = 2,colour="#000000")+
    geom_hline(yintercept = c(-log10(0.05)), lty = 2,colour="#000000")+
    xlab(NULL)+ylab(NULL)+
    theme_bw()+
    theme(
      panel.grid =element_blank(),
      plot.title = element_text(color = 'black',
                                hjust=0.5)
    )+annotate("text",x = c(-3, -3, 3, 3), y = c(0, -0.5, 0, -0.5), 
               label = c(paste0('All: ',alldown), paste0('RBP: ',RBPdown), 
                         paste0('All: ',allup), paste0('RBP: ',RBPup)))
  
  pdf(paste0('test/',compare0,'.pdf'))
  print(p)
  dev.off()
  try(dev.off())
  plotall[[i]]=p
}
library(patchwork)
pdf('FigureS1_A.pdf',width=10,height=8)
plotall[[1]]+plotall[[2]]+plotall[[3]]+plotall[[4]]+plotall[[5]]+plotall[[6]]+plot_layout(ncol=3,guides='collect')
dev.off()


#Protein

plotall=list()

for (i in 1:length(compare)) {
  compare0=compare[i]
  res2=read.csv(paste0('data/DEresult/Protein_',compare0,'_DE.txt'),sep='\t',header=T)
  res2$adj.P.Val=as.numeric(res2$adj.P.Val)
  res2$logFC=as.numeric(res2$logFC)
  threshold_new_res2 <- as.factor(ifelse(res2$adj.P.Val < 0.05 &
                                           abs(res2$logFC) > 0  ,
                                         ifelse(res2$logFC > 0 ,'Up','Down'),'Not'))
  res2$threshold=as.character(threshold_new_res2)
  alldown=length(which(res2$threshold == 'Down'))
  allup=length(which(res2$threshold == 'Up'))
  
  res2$Size=1
  res2[which(rownames(res2) %in% RBP & res2$threshold != 'Not'),'Size'] =2
  res2$alpha=0.5
  res2[which(rownames(res2) %in% RBP & res2$threshold != 'Not'),'alpha'] =0.8
  res2$shape='Other'
  res2[which(rownames(res2) %in% RBP & res2$threshold != 'Not'),'shape'] ='RBP'
  RBPdown=length(which(res2$threshold == 'Down' & res2$shape == 'RBP'))
  RBPup=length(which(res2$threshold == 'Up' & res2$shape == 'RBP'))
  
  
  p=ggplot(res2,aes(x=logFC,y=-log10(adj.P.Val),colour=threshold)) +
    xlab("log2(Fold Change)")+ylab("-log10(padjust)") +ggtitle(compare0)+
    geom_point(aes(alpha=alpha,size = Size,shape = shape))+
    scale_size_continuous(range = c(2,3))+
    
    xlim(-4,4) +
    scale_color_manual(values =c("#0072B5","grey","#BC3C28"))+
    geom_vline(xintercept = c(0), lty = 2,colour="#000000")+
    geom_hline(yintercept = c(-log10(0.05)), lty = 2,colour="#000000")+
    xlab(NULL)+ylab(NULL)+
    
    theme_bw()+
    theme(
      panel.grid =element_blank(),
      plot.title = element_text(color = 'black',
                                hjust=0.5)
    )+annotate("text",x = c(-3, -3, 3, 3), y = c(0, -0.5, 0, -0.5), 
               label = c(paste0('All: ',alldown), paste0('RBP: ',RBPdown), 
                         paste0('All: ',allup), paste0('RBP: ',RBPup)))
  pdf(paste0('test/',compare0,'.pdf'))
  print(p)
  dev.off()
  try(dev.off())
  plotall[[i]]=p
  
}
pdf('FigureS1_B.pdf',width=10,height=8)
plotall[[1]]+plotall[[2]]+plotall[[3]]+plotall[[4]]+plotall[[5]]+plotall[[6]]+plot_layout(ncol=3,guides='collect')
dev.off()