outloc = commandArgs(trailingOnly=TRUE)[1] 
options(stringsAsFactors = FALSE)
library(tidyverse)
library(ggpubr)
library(broom)


res=read_delim(paste0(outloc,"/binding_site_removal/coexpression.txt.gz"),
               delim = "\t",col_names = F,
               col_types = cols(
                 X1 = col_character(),
                 X2 = col_double(),
                 X3 = col_double(),
                 X4 = col_character(),
                 X5 = col_character()
               ))

lapply(res$X1%>%unique()%>%.[!grepl("no_mutation",.)], function(y){
  res=res%>%filter(X1%in%c(y,"no_mutation"))
  res$X4%>%str_split(.,pattern = ",")%>%unlist()%>%unique()%>%as.numeric()%>%
    sort()%>%length -> n_length
  
  RNA_design=matrix(0,nrow=nrow(res),ncol=n_length)
  colnames(RNA_design)=paste0("RNA_design",1:n_length)
  for(i in 1:nrow(res)){
    indx=res$X4[i]%>%str_split(.,pattern = ",")%>%unlist()%>%unique()%>%as.numeric() + 1
    RNA_design[i,indx]=1
  }
  
  res$X5%>%str_split(.,pattern = ",")%>%unlist()%>%unique()%>%as.numeric()%>%
    sort()%>%length -> n_length
  
  DNA_design=matrix(0,nrow=nrow(res),ncol=n_length)
  colnames(DNA_design)=paste0("DNA_design",1:n_length)
  for(i in 1:nrow(res)){
    indx=res$X5[i]%>%str_split(.,pattern = ",")%>%unlist()%>%unique()%>%as.numeric() + 1
    DNA_design[i,indx]=1
  }
  
  lm(X2~.,data =cbind(RNA_design,DNA_design)%>%
       as.data.frame()%>%
       mutate(X2=res$X2))%>%broom::tidy()%>%
    filter(term!="(Intercept)")%>%
    mutate(term=gsub("RNA_design","RNA_interval_",term),
           term=gsub("DNA_design","promoter_interval_",term))%>%
    mutate(Gene=y)%>%return()
})%>%bind_rows()%>%
  arrange(p.value) -> lm_res

lm_res%>%
  write.table(.,file = paste0(outloc,"/binding_site_removal/regression_res.txt"),
              append = F,quote = F,sep = "\t",row.names = F,col.names = T)

