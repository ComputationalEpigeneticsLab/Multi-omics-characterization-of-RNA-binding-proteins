
SARS3=read.csv('data/SARS_RBP_intersect.csv',header=T)
SARS3_RBP=unique(SARS3$Gene)
dat_td=read.csv('Drug/clear_all_TTD_drug_target.txt',sep='\t',header=T)
test=read.csv('Drug/final_Z_P_0.001_2.txt',sep='\t',header=T)
test=merge(test,dat_td,by.x='Drug',by.y='DrugID')
test_2=test[grep('Approved',test$DRUGCLAS),]
drug=test_2[order(test_2$Z_score),]
drug2=drug[,c('Drug','Z_score','GENENAME.x','DrugName','MOA')] %>% distinct()
colnames(drug2)[3]='Target'


humanP_humanP=read.csv('data/humanP_humanP.txt',sep='\t',header=T)
humanP_humanP_1=humanP_humanP %>% filter(InteractionA %in% drug2$Target)
humanP_humanP_2=humanP_humanP %>% filter(InteractionB %in% drug2$Target)
humanP_humanP_2=humanP_humanP_2[,c(2,1)]
humanP_humanP2=rbind(humanP_humanP_1,humanP_humanP_2)
colnames(humanP_humanP2)=c('Target','Symbol')

sarscov=read.csv('data/SARS-CoV-2.txt',sep='\t',header=T,check.names=F)
# write.table(sarscov$'Gene-2','data/SARS-CoV-2_uniprotid.txt',sep='\t',quote=F,col.names=F,row.names=F)
transid=read.csv('data/SARS-CoV-2_uniprotid_transsymbol.txt',sep='\t',header=T,check.names=F)
sarscov2=merge(sarscov,transid,by.x='Gene-2',by.y='From')
sarscov2=sarscov2[,c('Gene1','To')]
colnames(sarscov2)=c('Virus','Symbol')
extrainfo=openxlsx::read.xlsx('data/41422_2021_581_MOESM8_ESM.xlsx',sheet=2)
extrainfo$Gene1='SARS-CoV-2'
extrainfo=extrainfo[,c(ncol(extrainfo),2)]
colnames(extrainfo)=c('Virus','Symbol')
extrainfo=extrainfo[-which(extrainfo$Symbol %in% sarscov2$Symbol),]
sarscov_all=rbind(sarscov2,extrainfo)

drug2_test=drug2[order(drug2$Z_score),]
select_drug = unique(drug2_test$Drug)[1:50]
drug2_test=drug2_test[which(drug2_test$Drug %in% select_drug),]

drug2_test2=merge(drug2_test,humanP_humanP2,by='Target',all.x=TRUE)

drug2_test2=drug2_test2 %>% filter(Symbol %in% SARS3_RBP)
drug2_test2=drug2_test2 %>% dplyr::select(DrugName,Target,Symbol)
colnames(drug2_test2)=c('Drugs','Targets','SARS-Cov-2 related RBPs')

data=drug2_test2
data$flow<-1:nrow(data)
data2<-reshape2::melt(data,id='flow')
variable<-summary(data2$variable)
data2$link<-1

circos_0=data[,c('Drugs','Targets')] %>% distinct()
colnames(circos_0)=c('Node_1','Node_2')
circos_1=data[,c('SARS-Cov-2 related RBPs','Targets')] %>% distinct()
colnames(circos_1)=c('Node_1','Node_2')
circos_prepare=rbind(circos_0,circos_1)

NodeType=data2[,c('value','variable')] %>% distinct()
colnames(NodeType)=c('Node','Type')
NodeType$Node=as.character(NodeType$Node)
NodeType$Type=as.character(NodeType$Type)

brand = c(structure(NodeType$Type, names=NodeType$Node))
brand = brand[!duplicated(names(brand))]
brand = brand[order(brand, names(brand))]
brand=factor(brand,levels = c('Drugs','Targets','SARS-Cov-2 related RBPs'))
brand=sort(brand)


# color
library(RColorBrewer)
mycolor<-c('#96c08e','#75a7c4','#dd8994')
brand_color = structure(mycolor, names = levels(brand))
model_color = structure(rep(mycolor,as.numeric(table(brand))), names = names(brand))
pdf('Drug/FigureS7A-Drug_Target_RBP_TTD_circos2.pdf')
library(circlize)
circos.clear()
gap.after = do.call("c", lapply(table(brand), function(i) c(rep(2, i-1), 8)))
circos.par(gap.after = gap.after, cell.padding = c(0, 0, 0, 0))

chordDiagram(circos_prepare, order = names(brand), grid.col = model_color,transparency = 0.2,
             directional = 0, 
             annotationTrack = "grid", preAllocateTracks = list(
               list(track.height = 0.02)
               ),
             annotationTrackHeight = mm_h(c(4, 2))
)

circos.track(track.index = 2, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.index = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), mean(ylim), sector.index, col = "black", cex =0.4, 
              facing = "reverse.clockwise", niceFacing = T)
}, bg.border = NA)

for(b in unique(brand)) {
  model = names(brand[brand == b])
  highlight.sector(sector.index = model, track.index = 1, col = brand_color[b],
                   text = b, text.vjust = -1, niceFacing = TRUE
                   )
}

circos.clear()
dev.off()


SARS3=read.csv('data/SARS_RBP_intersect.csv',header=T)
SARS3_RBP=unique(SARS3$Gene)
drug=read.csv('Drug/final_Z_P_0.001_2.txt',sep='\t',header=T)
drug2=drug[,c(8,3,6)]
colnames(drug2)[3]='Target'

humanP_humanP=read.csv('data/humanP_humanP.txt',sep='\t',header=T)
humanP_humanP_1=humanP_humanP %>% filter(InteractionA %in% drug2$Target)
humanP_humanP_2=humanP_humanP %>% filter(InteractionB %in% drug2$Target)
humanP_humanP_2=humanP_humanP_2[,c(2,1)]
humanP_humanP2=rbind(humanP_humanP_1,humanP_humanP_2)
colnames(humanP_humanP2)=c('Target','Symbol')


#Doxorubicin
drug2_test=drug2[order(drug2$Z_score),]
drug2_test2=merge(drug2_test,humanP_humanP2,by='Target')

drug2_test3=merge(drug2_test2,SARS3,by.x='Symbol',by.y='Gene')
drug2_test3=drug2_test3 %>% dplyr::select(DrugName,Target,Symbol,Gene1,Z_score)
colnames(drug2_test3)=c('Drug','Target','SARS_RBP','SARS_Cov','Z_score')
drug2_test3=drug2_test3[order(drug2_test3$Z_score),]
drug0='Doxorubicin'
drug2_test4=drug2_test3 %>% filter(Drug ==drug0)
drug2_test4=drug2_test4[!duplicated(drug2_test4),]

drug_info=cbind(drug2_test4$Drug,'Drug')
target_info=cbind(drug2_test4$Target,'Target')
SARS_RBP_info=cbind(drug2_test4$SARS_RBP,'SARS_RBP')
SARS_Cov_info=cbind(drug2_test4$SARS_Cov,'SARS_Cov')
Nodetype=rbind(drug_info,target_info,SARS_RBP_info,SARS_Cov_info)
colnames(Nodetype)=c('Node','Type')
Nodetype=as.data.frame(Nodetype) %>% distinct()
write.table(Nodetype,paste0('Drug/NodeType_',drug0,'.txt'),sep='\t',quote=F,row.names=F,col.names=T)

drug2_test4_1=drug2_test4[,c(1,2)] %>% distinct() %>% as.matrix()
drug2_test4_1=cbind(drug2_test4_1,'drug')
drug2_test4_2=drug2_test4[,c(2,3)] %>% distinct() %>% as.matrix()
drug2_test4_2=cbind(drug2_test4_2,'target')
drug2_test4_3=drug2_test4[,c(3,4)] %>% distinct() %>% as.matrix()
drug2_test4_3=cbind(drug2_test4_3,'sars')
drug2_test5=rbind(drug2_test4_1,drug2_test4_2,drug2_test4_3)
colnames(drug2_test5)=c('Drug','Target','SARS_Cov')
write.table(drug2_test5,paste0('Drug/DrugNetwork_',drug0,'.txt'),sep='\t',quote=F,row.names=F,col.names=T)
#these file were submitted to Cytoscape for visualization

#Mycophenolate mofetil
drug2_test=drug2[order(drug2$Z_score),]
drug2_test2=merge(drug2_test,humanP_humanP2,by='Target')
drug2_test3=merge(drug2_test2,SARS3,by.x='Symbol',by.y='Gene')
drug2_test3=drug2_test3 %>% dplyr::select(DrugName,Target,Symbol,Gene1,Z_score)
colnames(drug2_test3)=c('Drug','Target','SARS_RBP','SARS_Cov','Z_score')
drug2_test3=drug2_test3[order(drug2_test3$Z_score),]
drug0='Mycophenolate mofetil'
drug2_test4=drug2_test3 %>% filter(Drug ==drug0)
drug2_test4=drug2_test4[!duplicated(drug2_test4),]

drug_info=cbind(drug2_test4$Drug,'Drug')
target_info=cbind(drug2_test4$Target,'Target')
SARS_RBP_info=cbind(drug2_test4$SARS_RBP,'SARS_RBP')
SARS_Cov_info=cbind(drug2_test4$SARS_Cov,'SARS_Cov')
Nodetype=rbind(drug_info,target_info,SARS_RBP_info,SARS_Cov_info)
colnames(Nodetype)=c('Node','Type')
Nodetype=as.data.frame(Nodetype) %>% distinct()
write.table(Nodetype,paste0('Drug/NodeType_',drug0,'.txt'),sep='\t',quote=F,row.names=F,col.names=T)

drug2_test4_1=drug2_test4[,c(1,2)] %>% distinct() %>% as.matrix()
drug2_test4_1=cbind(drug2_test4_1,'drug')
drug2_test4_2=drug2_test4[,c(2,3)] %>% distinct() %>% as.matrix()
drug2_test4_2=cbind(drug2_test4_2,'target')
drug2_test4_3=drug2_test4[,c(3,4)] %>% distinct() %>% as.matrix()
drug2_test4_3=cbind(drug2_test4_3,'sars')
drug2_test5=rbind(drug2_test4_1,drug2_test4_2,drug2_test4_3)
colnames(drug2_test5)=c('Drug','Target','SARS_Cov')
write.table(drug2_test5,paste0('Drug/DrugNetwork_',drug0,'.txt'),sep='\t',quote=F,row.names=F,col.names=T)
#these file were submitted to Cytoscape for visualization
