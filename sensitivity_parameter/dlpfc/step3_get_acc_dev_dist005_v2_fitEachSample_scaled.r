
### load pacakges 
library(tidyverse)
library(ggplot2)
library(Matrix)
library(ggforce)#
library(cowplot)
library(RColorBrewer)
library(grid)
require(nnet)

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

# setwd("/dcs04/hansen/data/ywang/ST/DLPFC/PASTE_out_dist005_scMean/")
# dir_processedData = "/dcs04/hansen/data/ywang/ST/DLPFC/processed_Data_dist005/"

################# 
dir_output="/dcs04/hansen/data/ywang/ST/disp/"
dir.create(dir_output)
setwd(dir_output)
#######################################
################# paste results
#######################################
list_filenames = c('0001', '001', '01', 1, 10, 100)
list_disp = c(0.001, 0.01, 0.1, 1, 10, 100)

list_acc_c = list()
list_dev_c = list()
for(k in c(1:6)){
  disp= list_disp[k]
  disp_name = list_filenames[k]
  
  for(ksample in 1:12){
  
    # load layers
    load(paste0("//dcs04/hansen/data/ywang/archive/scRNA/Oct5_2021_Lukas_data_more_Genes/out/layer_sample_Sample",ksample,".RData"))
    layer
    
    # load mNSF factors
    a_=read.csv(paste0("//dcs04/hansen/data/ywang/ST/DLPFC/PASTE_out_keepAll_scTransform/factors_nb_szMean_sample_s",
                         ksample,"_L10_fullData_disp",disp_name,".csv"),header=T)

    
    a = scale(a_[,-1])
    df_ = data.frame(a)
    df_$layer = layer
    
    df_=df_[!is.na(df_$layer),]
    
  
    multinom_model <- multinom(layer ~ ., data = df_)
    dev1=mean(multinom_model$deviance)
  
    # df_=df_[!is.na(df_$layer),]
    p.fit <- predict(multinom_model, predictors=df_, type='probs')
    pred_=unlist(lapply(1:nrow(p.fit),function(xx){
      colnames(p.fit)[which.max(p.fit[xx,])]
    }))
    # 
    acc1=mean(pred_==df_$layer,na.rm=T)
    
    if(ksample == 1){
      acc_c = acc1
      dev_c = dev1
    }else{
      acc_c = c(acc_c, acc1)
      dev_c = c(dev_c, dev1)
    }
    # 
  }
  list_acc_c[[as.character(disp_name)]] = acc_c
  list_dev_c[[as.character(disp_name)]] = dev_c
}

for(k in c(1:6)){
  disp= list_disp[k]
  disp_name = list_filenames[k]
  # print(nfactor)
  print(mean(list_acc_c[[as.character(disp_name)]]))
  print(mean(list_dev_c[[as.character(disp_name)]]))
}
saveRDS(list_acc_c, file = "list_acc_c_scaled.rds")
saveRDS(list_dev_c, file = "list_dev_c_scaled.rds")
