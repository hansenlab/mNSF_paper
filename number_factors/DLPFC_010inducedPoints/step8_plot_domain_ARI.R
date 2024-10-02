# ml conda_R/4.4
library(SpatialPCA)
library(tidyverse)
library(ggplot2)
library(Matrix)
library(ggforce)#
library(cowplot)
library(RColorBrewer)
library(grid)
require(nnet)
library(aricode)
library(matrixStats)


#######
dir_output="/dcs04/hansen/data/ywang/ST/DLPFC/threeSamples/"
setwd(dir_output)
dir_processedData = "/dcs04/hansen/data/ywang/ST/DLPFC/processed_Data_keepAll/"

#######################################
location_list = list()
for(ksample in 1:12){
  X=read.csv(paste0('//dcs04/hansen/data/ywang/ST/DLPFC/processed_Data///X_allSpots_sample',ksample,'.csv'))
  location_list[[ksample]] = X
}
list_layer=list()
for(ksample in 1:12){
  
  load(paste0("//dcs04/hansen/data/ywang/archive/scRNA/Oct5_2021_Lukas_data_more_Genes/out/layer_sample_Sample",ksample,".RData"))
  
  list_layer[[ksample]]=layer
  
}





# L=4,6,8,10
L=2
L=4
L=6
L=8
L=10

#######################################
for(kdonor in c(1:3)){
  for(k in 1:4){
    ksample = (kdonor-1)*4+ k
  }  
  if(L!=10){
    a_=read.csv(paste0("//dcs04/hansen/data/ywang/ST/DLPFC/PASTE_out_keepAll_scTransform//factors_nb_szMean_sample_s",
                       ksample,"_L",L,"_fullData_disp1.csv"),header=T)
  }else{
    a_=read.csv(paste0("//dcs04/hansen/data/ywang/ST/DLPFC/PASTE_out_keepAll_scTransform//factors_nb_szMean_sample_s",
                       ksample,"_L",L,"_fullData_disp001.csv"),header=T)
  }
  a = a_[,-1]
  if(kdonor ==1){a_c =a_}else{
    a_c = rbind(a_c,a_)
  }

}

sd_cols = colSds(data.matrix(a_c))


##################
for(kdonor in c(1:3)){
  for(k in 1:4){
    ksample = (kdonor-1)*4+ k
    layers=list_layer[[ksample]]
    
    ##
    
    if(L!=10){
      a_=read.csv(paste0("//dcs04/hansen/data/ywang/ST/DLPFC/PASTE_out_keepAll_scTransform//factors_nb_szMean_sample_s",
                         ksample,"_L",L,"_fullData_disp1.csv"),header=T)
    }else{
      a_=read.csv(paste0("//dcs04/hansen/data/ywang/ST/DLPFC/PASTE_out_keepAll_scTransform//factors_nb_szMean_sample_s",
                         ksample,"_L",L,"_fullData_disp001.csv"),header=T)
    }
    a = a_[,-1]
    X=read.csv(paste0('//dcs04/hansen/data/ywang/ST/DLPFC/processed_Data///X_allSpots_sample',ksample,'.csv'))
    
    
    layers_pred_sampleTmp = walktrap_clustering(clusternum=ifelse(kdonor == 2, 5, 7), latent_dat=t((a_[,]))/sd_cols, knearest=70)
    which_na = which(is.na(layers))
    print(ARI(layers_pred_sampleTmp[-which_na], layers[-which_na]))
    
    
    saveRDS(layers_pred_sampleTmp, file = paste0("layers_pred_withinEachsample_sample_",ksample,"_L",L,"_scaleV2.rds"))
    
  }
 
  
}


#######################################
ARI_vec=c()

for(kdonor in c(1:3)){
  for(k in 1:4){
    ksample = (kdonor-1)*4+ k
    layers=list_layer[[ksample]]
    # layers=list_layer[[ksample]]
  
  ##

  X=read.csv(paste0('//dcs04/hansen/data/ywang/ST/DLPFC/processed_Data///X_allSpots_sample',ksample,'.csv'))


  layers_pred_sampleTmp = readRDS(paste0("layers_pred_withinEachsample_sample_",ksample,"_L",L,"_scaleV2.rds"))

  clusterlabel_refine = refine_cluster_10x(clusterlabels=layers_pred_sampleTmp,
                                           location=location_list[[ksample]],shape="hexagon")
  
  which_na = which(is.na(layers))

  ART_tmp = ARI(clusterlabel_refine[-which_na], layers[-which_na])
  print(ART_tmp)
  ARI_vec = c(ARI_vec, ART_tmp)
  
  saveRDS(clusterlabel_refine, file = paste0("layers_pred_withinEachsample_sample_",ksample,"_L",L,"_refined_scaleV2.rds"))
  
}
}
saveRDS(ARI_vec,file = paste0("ARI_vec_sample1_5_9_predLayers_mnsf_ClusterWithinEachSample_mnsfrefined__L",L,"_scaleV2.rds"))
