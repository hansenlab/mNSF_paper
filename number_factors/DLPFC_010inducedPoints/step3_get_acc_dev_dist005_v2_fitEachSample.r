
### load pacakges 
library(tidyverse)
library(ggplot2)
library(Matrix)
library(ggforce)#
library(cowplot)
library(RColorBrewer)
library(grid)
require(nnet)


################# 
dir_output="/dcs04/hansen/data/ywang/ST/nfactor/"
setwd(dir_output)
#######################################
#######################################
list_factor = c(4,8,12,16,20)
list_acc_c = list()
list_dev_c = list()
for(nfactor in list_factor){
  for(ksample in 1:12){
  
    # load layers
    load(paste0("//dcs04/hansen/data/ywang/archive/scRNA/Oct5_2021_Lukas_data_more_Genes/out/layer_sample_Sample",ksample,".RData"))
    # layer

    a_=read.csv(paste0("//dcs04/hansen/data/ywang/ST/DLPFC/PASTE_out_keepAll_scTransform//factors_nb_szMean_sample_s",
                         ksample,"_L",nfactor,"_fullData_010Induced.csv"),header=T)

    
    a = a_[,-1]*10
    df_ = data.frame(a)
    df_$layer = layer

    multinom_model <- multinom(layer ~ ., data = df_)
    dev1=mean(multinom_model$deviance)
  
    df_=df_[!is.na(df_$layer),]
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
  list_acc_c[[as.character(nfactor)]] = acc_c
  list_dev_c[[as.character(nfactor)]] = dev_c
  
  # print(nfactpr)
  
  # print(mean(acc_c))
  # print(mean(dev_c))
  
}
for(nfactor in list_factor){
  print(nfactor)
  print(mean(list_acc_c[[as.character(nfactor)]]))
  print(mean(list_dev_c[[as.character(nfactor)]]))
  
}
saveRDS(list_acc_c, file = "list_acc_c_010Induced.rds")
saveRDS(list_dev_c, file = "list_dev_c_010Induced.rds")

