
# Descriptions of outputs
# RNA_importance_mean.txt: Average DeepLIFT score for each RNA regulator in each tissue.
# promoter_importance_mean.txt: Average DeepLIFT score for each promoter regulator in each tissue.

outloc = commandArgs(trailingOnly=TRUE)[1] 
options(stringsAsFactors = FALSE)
library(tidyverse)
library(ggpubr)
library(broom)


feature_names = read_delim(paste0(dirname(outloc),"/feature_norm_stats.txt"),
                           delim = ",",col_names = T)

# TF importance ####
feature_names%>%
  filter(feature_type=="promoter_range")%>%
  .$feature_name%>%
  gsub("promoter_annot_|_GTRD.rds|.rds","",.)-> col_name
feature_names%>%
  filter(feature_type=="deg_stat")%>%
  arrange(row_indx) -> deg_stat

lapply(0:(nrow(deg_stat)-1), function(sample_indx){
  read_delim(paste0(outloc,"/DeepLIFT/DNA_",sample_indx,".txt.gz"),
             delim = "\t",col_names = F) -> dLIFT_ds
  
  dLIFT_mat=matrix(0,nrow=nrow(dLIFT_ds),ncol = length(col_name))
  rownames(dLIFT_mat)=dLIFT_ds$X2
  colnames(dLIFT_mat)=col_name
  
  for (i in 1:nrow(dLIFT_ds)){
    dLIFT_mat[i,]=dLIFT_ds%>%.$X3%>%.[i]%>%
      str_split(.,pattern = ",")%>%.[[1]]%>%as.numeric()
  }
  
  dLIFT_mat%>%as.data.frame()%>%
    mutate(Gene=rownames(dLIFT_mat))%>%
    gather(feat_name,DeepLIFT,-Gene)%>%
    mutate(sample_name=deg_stat$feature_name[sample_indx+1])%>%
    ungroup%>%
    group_by(feat_name,sample_name)%>%
    summarise(DeepLIFT_mean=mean(DeepLIFT))%>%
    dplyr::select(sample_name,feat_name,DeepLIFT_mean)%>%
    return()
})%>%bind_rows()%>%
  spread(feat_name,DeepLIFT_mean)->promoter_importance_mean
write.table(promoter_importance_mean,file = paste0(outloc,"/DeepLIFT/promoter_importance_mean.txt"),
            append = F,quote = F,sep = "\t",row.names = F,col.names = T)


# RNA importance ####
feature_names%>%
  filter(feature_type=="mRNA_range")%>%
  .$feature_name-> col_name
feature_names%>%
  filter(feature_type=="deg_stat")%>%
  arrange(row_indx) -> deg_stat

lapply(0:(nrow(deg_stat)-1), function(sample_indx){
  read_delim(paste0(outloc,"/DeepLIFT/RNA_",sample_indx,".txt.gz"),
             delim = "\t",col_names = F) -> dLIFT_ds
  
  dLIFT_mat=matrix(0,nrow=nrow(dLIFT_ds),ncol = length(col_name))
  rownames(dLIFT_mat)=dLIFT_ds$X2
  colnames(dLIFT_mat)=col_name
  
  for (i in 1:nrow(dLIFT_ds)){
    dLIFT_mat[i,]=dLIFT_ds%>%.$X3%>%.[i]%>%
      str_split(.,pattern = ",")%>%.[[1]]%>%as.numeric()
  }
  
  dLIFT_mat%>%as.data.frame()%>%
    mutate(Gene=rownames(dLIFT_mat))%>%
    gather(feat_name,DeepLIFT,-Gene)%>%
    mutate(sample_name=deg_stat$feature_name[sample_indx+1])%>%
    ungroup%>%
    group_by(feat_name,sample_name)%>%
    summarise(DeepLIFT_mean=mean(DeepLIFT))%>%
    dplyr::select(sample_name,feat_name,DeepLIFT_mean)%>%
    return()
})%>%bind_rows()%>%
  spread(feat_name,DeepLIFT_mean)->RNA_importance_mean
write.table(RNA_importance_mean,file = paste0(outloc,"/DeepLIFT/RNA_importance_mean.txt"),
            append = F,quote = F,sep = "\t",row.names = F,col.names = T)