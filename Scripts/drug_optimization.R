
#random network
rm(list = ls())
setwd("SARS_drug")
library(dplyr)
library(igraph)
humanP_humanP<-read.table("humanP_humanP.txt",stringsAsFactors = F,header = T)
humanP_humanP<-humanP_humanP %>% distinct()
g<-make_graph(t(humanP_humanP[,1:2]),directed = F)
shortest_path<-shortest.paths(g)

drug_target<-read.table("drug_humanP.txt",stringsAsFactors = F,header = T)
drug_target<-drug_target[order(drug_target$DrugID),]
drugID<-read.table("drugID.txt",stringsAsFactors = F,header = T)
drugID<-unique(drugID)


subnetwork2<-graph.edgelist(as.matrix(humanP_humanP[,1:2]),directed = F)
realnetwork<-get.adjacency(subnetwork2)
realnetwork<- as.matrix(realnetwork)
row.degree<- rowSums(realnetwork)
col.degree<- colSums(realnetwork)
row.degree<- as.matrix(row.degree)


cbind.all <- function (...) 
{
  nm <- list(...)
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow))
  do.call(cbind, lapply(nm, function(x) rbind(x, matrix(, n - 
                                                          nrow(x), ncol(x)))))
}
setwd("SARS_drug/result")

filename_select<-list.files("SARS_drug/virus_mcode_select")
Sys.time()

distance_random_all<-data.frame()

for(o in 1:1000){
  new=degree.sequence.game(row.degree, in.deg = NULL,method = "vl")  
  shortest_path1<-shortest.paths(new)
  rownames(shortest_path1)<-rownames(shortest_path)
  colnames(shortest_path1)<-colnames(shortest_path)
  
  virus_type_all<-c()
  drugID_all<-c()
  
  
  print(paste('o =',o))
  distance_random<-c()
  for (l in 1:length(filename_select)){
    select_mcode<-read.table(paste0("SARS_drug/virus_mcode_select/",
                                    filename_select[l]),as.is=T,header = T,sep = "\t")
    print(paste('l =',l))
    
    
    net<-c()
    for(m in 1:length(drugID$DrugID)){
      drugtarget<-drug_target[which(drug_target$DrugID %in% drugID$DrugID[m]),2]
      
      virus_type<-gsub('_(.*)','',filename_select[l])
      virus_type_all<-c(virus_type_all,virus_type)
      drugID_all<-c(drugID_all,drugID$DrugID[m])
      
      
      virus_sum=0
      virushumanP_num<-length(select_mcode$symbol)
      virus_min=shortest_path1[which(rownames(shortest_path1) %in% select_mcode[,1]),
                               which(colnames(shortest_path1) %in% drugtarget)]
      if(class(virus_min)[1]=="numeric"){
        virus_sum=virus_sum+sum(virus_min)
      }
      if(class(virus_min)[1]!="numeric"){
        virus_min<-apply(virus_min,2,min)
        virus_sum=virus_sum+sum(virus_min)
      }
      
      
      drug_sum=0    
      drugtarget_num<-length(drugtarget)
      drug_min=shortest_path1[which(rownames(shortest_path1) %in% select_mcode[,1]),
                              which(colnames(shortest_path1) %in% drugtarget)]
      if(class(drug_min)[1]=="numeric"){
        drug_sum=drug_sum+min(drug_min)
      }
      if(class(drug_min)[1]!="numeric"){
        drug_min<-apply(drug_min,1,min)
        drug_sum=drug_sum+sum(drug_min)
      }
      
      network_proximity<-(virus_sum+drug_sum)/(virushumanP_num+drugtarget_num)
      distance_random<-c(distance_random,network_proximity)
      
    }
    
    
    
  }
  distance_random_all<-cbind.all(distance_random_all,distance_random)
  distance_random_all<-as.data.frame(distance_random_all)
  cat(o,sep = "\n")
  
  
  
}


distance_random_all2<-cbind(cbind(virus_type_all,drugID_all),distance_random_all)
colnames(distance_random_all2)[3:(3+o-1)]<-rep(paste0('Number_',1:o))
write.table(distance_random_all2,'distance_random_all_2_1.txt',sep='\t',quote=F,row.names = F)


print('End!!')


##############################################################
#Real network
options(repos = structure(c(CRAN="https://mirror.tuna.tsinghua.edu.cn/CRAN/")))
setwd("SARS_drug/")
library(dplyr)
library(igraph)
humanP_humanP<-read.table("humanP_humanP.txt",stringsAsFactors = F,header = T)
humanP_humanP<-humanP_humanP %>% distinct()

g<-make_graph(t(humanP_humanP[,1:2]),directed = F)
shortest_path<-shortest.paths(g)

drug_target<-read.table("drug_humanP.txt",stringsAsFactors = F,header = T)
drug_target<-drug_target[order(drug_target$DrugID),]
drugID<-read.table("drugID.txt",stringsAsFactors = F,header = T)
drugID<-unique(drugID)


filename_select <- list.files("SARS_drug/virus_mcode_select")
setwd("SARS_drug/distance_final")

file.create("distance.txt")

Sys.time()
for (l in 1:length(filename_select)){
  select_mcode<-read.table(paste0("SARS_drug/virus_mcode_select/",
                                  filename_select[l]),as.is=T,header = T,sep = "\t")
  
  net<-c()
  Sys.time()
  for(m in 1:length(drugID$DrugID)){
    drugtarget<-drug_target[which(drug_target$DrugID %in% drugID$DrugID[m]),2]
    
    
    virus_sum=0
    virushumanP_num<-length(select_mcode$symbol)
    virus_min=shortest_path[which(rownames(shortest_path) %in% select_mcode[,1]),
                            which(colnames(shortest_path) %in% drugtarget)]
    if(class(virus_min)[1]=="numeric"){
      virus_sum=virus_sum+sum(virus_min)
    }
    if(class(virus_min)[1]!="numeric"){
      virus_min<-apply(virus_min,2,min)
      virus_sum=virus_sum+sum(virus_min)
    }
    
    
    drug_sum=0    
    drugtarget_num<-length(drugtarget)
    drug_min=shortest_path[which(rownames(shortest_path) %in% select_mcode[,1]),
                           which(colnames(shortest_path) %in% drugtarget)]
    if(class(drug_min)[1]=="numeric"){
      drug_sum=drug_sum+min(drug_min)
    }
    if(class(drug_min)[1]!="numeric"){
      drug_min<-apply(drug_min,1,min)
      drug_sum=drug_sum+sum(drug_min)
    }
    
    network_proximity<-(virus_sum+drug_sum)/(virushumanP_num+drugtarget_num)
    net<-c(net,network_proximity)
    
    myfile_all<-gsub('_(.*)','',filename_select[l])
    distance<-t(c(myfile_all,drugID$DrugID[m],network_proximity))
    write.table(distance,"distance.txt",row.names = F,
                col.names = F, sep = "\t", append=T, quote = F)
    
    
  }
  
  
  cat(l,sep = "\n")
  
  
}
Sys.time()




###############################################################################

rm(list = ls())
###random network
setwd("SARS_drug/result")
distance_random_all_1<-read.table("distance_random_all_2_1.txt",stringsAsFactors = F,header = T)
distance_random_all_2<-read.table("distance_random_all_2_2.txt",stringsAsFactors = F,header = T)
distance_random_all_3<-read.table("distance_random_all_2_3.txt",stringsAsFactors = F,header = T)
distance_random_all_4<-read.table("distance_random_all_2_4.txt",stringsAsFactors = F,header = T)

distance_random_all_2<-distance_random_all_2[,-c(1:2)]
distance_random_all_3<-distance_random_all_3[,-c(1:2)]
distance_random_all_4<-distance_random_all_4[,-c(1:2)]
distance_random_final<-cbind(distance_random_all_1,distance_random_all_2,
                             distance_random_all_3,distance_random_all_4)
colnames(distance_random_final)[3:1002]<-rep(paste0('Number_',1:1000))

#####mean, sd
distance_random_final$mean<-apply(distance_random_final[3:1002],1,mean)
distance_random_final$standard<-apply(distance_random_final[3:1002],1,sd)
distance_random_final<-distance_random_final[-which(distance_random_final$Number_1=="Inf"),]

###Real network
setwd("SARS_drug/distance_final")
distance_actual<-read.table("distance.txt",stringsAsFactors = F,header = F)
colnames(distance_actual)<-c("virus_type_all","drugID_all","distance_value")
distance_actual<-distance_actual[-which(distance_actual$distance_value=="Inf"),]

file.create('final_value_big.txt')

setwd("SARS_drug/distance_final")
Sys.time()
####Z_score, P
for (k in 1:length(distance_actual[,1])) {
  
  ####z_score
  distance_exact<-distance_actual$distance_value
  distance_mean<-distance_random_final$mean
  distance_std<-distance_random_final$standard
  virus<-distance_actual$virus_type_all
  drug<-distance_actual$drugID_all
  Z_score<-(distance_exact-distance_mean)/distance_std
  
  
  less_num<-length(which(distance_random_final[k,3:1002]<distance_actual[k,3]))
  
  P_value<-less_num/1000
  
  final_value<-t(c(virus[k],drug[k],Z_score[k],P_value))
  write.table(final_value,"final_value_big.txt",row.names = F,
              col.names = F, sep = "\t", append=T, quote = F)
  
  
}
Sys.time()

setwd("SARS_drug/distance_final")
final_value_out<-read.table("final_value_big.txt",stringsAsFactors = F,header = F)
colnames(final_value_out)<-c("virus","Drug","Z_score","P_value")

p<-p.adjust(final_value_out[,4],n=nrow(final_value_out),method = "bonferroni")
index1<-which(final_value_out[,4]<0.001)
index2<-which(p<0.001)
final_value_out_adjust<-final_value_out[index2,]
final_value_out_adjust2<-final_value_out_adjust[which(final_value_out_adjust$Z_score< -1.5),]

write.table(final_value_out_adjust2,"final_Z_P_0.001.txt",row.names=F,col.names=T,quote=F,sep="\t")










