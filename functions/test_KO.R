outloc = commandArgs(trailingOnly=TRUE)[1] 
options(stringsAsFactors = FALSE)
library(tidyverse)
library(ggpubr)
library(broom)


feature_names = read_delim(paste0(dirname(outloc),"/feature_norm_stats.txt"),
                           delim = ",",col_names = T)

res=read_delim(paste0(outloc,"/regulator_KO/coexpression.txt.gz"),
               delim = "\t",col_names = F,col_types = cols(
                 X1 = col_double(),
                 X2 = col_double(),
                 X3 = col_character(),
                 X4 = col_character()
               ))

res$X3%>%str_split(.,pattern = ",")%>%unlist()%>%unique()%>%as.numeric()%>%
  sort()->feature_indx 

RNA_design=NULL
if (length(feature_indx)>0){
  feature_indx = data.frame(feature_indx)
  n_length=nrow(feature_indx)
  RNA_design=matrix(0,nrow=nrow(res),ncol=n_length)
  colnames(RNA_design)=feature_names%>%filter(feature_type=="mRNA_range")%>%
    filter(row_indx%in%feature_indx$feature_indx)%>%.$feature_name
  for(i in 1:nrow(res)){
    indx=feature_indx$feature_indx%in%(res$X3[i]%>%str_split(.,pattern = ",")%>%unlist()%>%unique()%>%as.numeric())
    RNA_design[i,indx]=1
  }
}

res$X4%>%str_split(.,pattern = ",")%>%unlist()%>%unique()%>%as.numeric()%>%
  sort()->feature_indx 

DNA_design=NULL
if (length(feature_indx)>0){
  feature_indx = data.frame(feature_indx)
  n_length=nrow(feature_indx)
  
  DNA_design=matrix(0,nrow=nrow(res),ncol=n_length)
  colnames(DNA_design)=feature_names%>%filter(feature_type=="promoter_range")%>%
    filter(row_indx%in%feature_indx$feature_indx)%>%.$feature_name%>%
    gsub("promoter_annot_|_GTRD.rds|.rds","",.)
  for(i in 1:nrow(res)){
    indx=feature_indx$feature_indx%in%(res$X4[i]%>%str_split(.,pattern = ",")%>%unlist()%>%unique()%>%as.numeric())
    DNA_design[i,indx]=1
  }
}

lm(res$X1~cbind(RNA_design,DNA_design))%>%broom::tidy()%>%
  mutate(term=gsub("cbind\\(RNA_design, DNA_design\\)","",term))%>%
  filter(term!="(Intercept)")%>%
  arrange(p.value)%>%
  write.table(.,file = paste0(outloc,"/regulator_KO/regression_res.txt"),
              append = F,quote = F,sep = "\t",row.names = F,col.names = T)
