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

library(clusterProfiler)
library(enrichplot)
genelist=unique(SARS3$Gene)
genelist=bitr(genelist,fromType='SYMBOL',toType=c('ENTREZID'),OrgDb='org.Hs.eg.db')
BP <- enrichGO(genelist$ENTREZID,"org.Hs.eg.db",ont="BP",
               keyType = "ENTREZID",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.1,
               readable = T)
genelist=unique(SARS3$Gene)
genelist=bitr(genelist,fromType='SYMBOL',toType=c('ENTREZID'),OrgDb='org.Hs.eg.db')
BP <- enrichGO(genelist$ENTREZID,"org.Hs.eg.db",ont="BP",
               keyType = "ENTREZID",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.1,
               readable = T)
pdf('GO_similarity_VirualRBP.pdf',width=14,height=10)
cnetplot(BP,showCategory = 10)
cnetplot(BP,showCategory = 10,circular=TRUE,colorEdge=TRUE)
barplot(BP, x = "GeneRatio", color = "p.adjust", showCategory = 10, size = NULL, split = NULL, font.size = 12, title="Dotplot for Gene Ontology Analysis")
goplot(BP)
heatplot(BP)
mat = GO_similarity(unique(BP$ID),ont='BP')
df = simplifyGO(mat,fontsize=c(10,20))
dev.off()
saveRDS(BP,'VirualRBP_BP.Rds')
genelist=setdiff(RBP,unique(SARS3$Gene))
genelist=bitr(genelist,fromType='SYMBOL',toType=c('ENTREZID'),OrgDb='org.Hs.eg.db')
BP <- enrichGO(genelist$ENTREZID,"org.Hs.eg.db",ont="BP",
               keyType = "ENTREZID",pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.1,
               readable = T)
pdf('GO_similarity_nonVirualRBP.pdf',width=14,height=10)
mat = GO_similarity(unique(BP$ID),ont='BP')
df = simplifyGO(mat,fontsize=c(10,20))
dev.off()