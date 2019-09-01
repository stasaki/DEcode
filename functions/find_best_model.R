
# Descriptions of outputs
#loss_parameters.rds: Hyper paramters and performance of the 10 training runs.
#best_model.txt: The best model among the 10 training runs. 

outdir = commandArgs(trailingOnly=TRUE)[1] 
#outdir = "./01_GTEX/network01/output/"
summary_outdir = paste0(outdir,"/summary/")
dir.create(summary_outdir)
options(stringsAsFactors = F)
library(tidyverse)
library(rjson)

# training history
list.files(path =outdir,pattern = "history.json",full.names = T)%>%
  lapply(., function(x){
    #print(x)
    ln_history=fromJSON(file = x)%>%
      as.data.frame()%>%
      mutate(run_time_stamp=basename(x)%>%gsub("_history.json","",.))
    
    ln_history=ln_history[,sort(colnames(ln_history))]
    return(ln_history)
  })%>%
  bind_rows()->ln_history

# paramters
list.files(path = outdir,pattern = "params.json",full.names = T)%>%
  lapply(., function(x){
    ln_params=fromJSON(file = x)%>%
      as.data.frame()%>%
      mutate(run_time_stamp=basename(x)%>%gsub("_params.json","",.))
    return(ln_params)
  })%>%
  bind_rows()->ln_params

# testing performance 
list.files(path = outdir,pattern = "*test_performance.txt",full.names = T)%>%
  lapply(., function(x){
    
    data.frame(run_time_stamp=basename(x)%>%gsub("_test.+","",.),
               test_loss = read.delim(x,header = F)%>%.$V1%>%.[1],
               test_pcor = read.delim(x,header = F)%>%.$V1%>%.[2])%>%return()
    
  })%>%bind_rows() -> test_performance

# save best loss for each parameter
ln_history%>%
  group_by(run_time_stamp)%>%
  arrange(val_loss)%>%
  filter(!duplicated(run_time_stamp))%>%
  dplyr::select(val_loss,val_pcor,run_time_stamp)%>%
  inner_join(.,ln_params,by="run_time_stamp")%>%
  inner_join(.,test_performance,by="run_time_stamp")->ds
saveRDS(ds,paste0(summary_outdir,"/loss_parameters.rds"))

#print("best model")
#print(ds[1,c(1:2,ncol(ds)-1,ncol(ds))])
#print(ds[1,3])

write.table(ds[1,3],file = paste0(summary_outdir,"/best_model.txt"),
            append = F,quote = F,sep = "\t",row.names = F,col.names = F)
