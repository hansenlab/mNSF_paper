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





set=2
set=3

set=4


#######################################

# Factors_df.to_csv(path.join(dir_output,"factors_nb_szMean_sample_s"+str(k+1)+"_L10_3samples_0.35Induced_svgs_set4.csv"))

##################
for(kdonor in c(1:3)){
  # for(k in 1){
    ksample = (kdonor-1)*4+ set 
    layers=list_layer[[ksample]]
    print(ksample)
    ##
    
    a_=read.csv(paste0("//dcs04/hansen/data/ywang/ST/DLPFC/threeSamples///factors_nb_szMean_sample_s",
                       kdonor,"_L10_3samples_0.35Induced_svgs_set",set,".csv"),header=T)
    a_ = a_[,-1]
    X=read.csv(paste0('//dcs04/hansen/data/ywang/ST/DLPFC/processed_Data///X_allSpots_sample',ksample,'.csv'))
    
    
    layers_pred_sampleTmp = walktrap_clustering(clusternum=ifelse(kdonor == 2, 5, 7), latent_dat=t(scale(a_[,])), knearest=70)
    which_na = which(is.na(layers))
    print(ARI(layers_pred_sampleTmp[-which_na], layers[-which_na]))
    
    
    saveRDS(layers_pred_sampleTmp, file = paste0("layers_pred_withinEachsample_sample_",ksample,"_scaleV1_3samples_010induced_set",set,".rds"))
    
  # }
  
  
}


#######################################
ARI_vec=c()

for(kdonor in c(1:3)){
  # for(k in 1){
    ksample = (kdonor-1)*4+ set 
    layers=list_layer[[ksample]]
    
    ##
    
    X=read.csv(paste0('//dcs04/hansen/data/ywang/ST/DLPFC/processed_Data///X_allSpots_sample',ksample,'.csv'))
    
    
    layers_pred_sampleTmp = readRDS(paste0("layers_pred_withinEachsample_sample_",ksample,"_scaleV1_3samples_010induced_set",set,".rds"))
    # saveRDS(layers_pred_sampleTmp, file = paste0("layers_pred_withinEachsample_sample_",ksample,"_scaleV1_3samples_010induced.rds"))
    
    clusterlabel_refine = refine_cluster_10x(clusterlabels=layers_pred_sampleTmp,
                                             location=location_list[[ksample]],shape="hexagon")
    
    which_na = which(is.na(layers))
    
    ART_tmp = ARI(clusterlabel_refine[-which_na], layers[-which_na])
    print(ART_tmp)
    ARI_vec = c(ARI_vec, ART_tmp)
    
    saveRDS(clusterlabel_refine, file = paste0("layers_pred_withinEachsample_sample_",ksample,"_refined_scaleV1_3samples_010induced_set",set,".rds"))
    
  # }
}
saveRDS(ARI_vec, file = paste0("ARI_vec_sample1_5_9_predLayers_mnsf_ClusterWithinEachSample_mnsfrefined_scaleV1_3samples_034induced_set",set,".rds"))
# set2
# [1] 0.3213066
# [1] 0.2665737
# [1] 0.3236757

# set3
# 1]  0.3455376
# [1] 0.5899907
# [1] 0.2348437

# set4
# [1] 0.2706302
# [1] 0.4954437
# [1] 0.2586318





