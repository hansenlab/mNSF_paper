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



# set=1

# set=2
# set=3

# set=4


#######################################
for(kdonor in c(1:3)){
  # for(k in 1){
  ksample = (kdonor-1)*4+ set 
  # }  
  if(set==1){    
    a_=read.csv(paste0("/dcs04/hansen/data/ywang/ST/DLPFC/threeSamples/factors_nb_szMean_sample_s",
                       kdonor,"_L10_3samples_0.35Induced_svgs.csv"),header=T)#factors for 3-sample mNSF
  }else{
    a_=read.csv(paste0("//dcs04/hansen/data/ywang/ST/DLPFC/threeSamples/factors_nb_szMean_sample_s",
                       kdonor,"_L10_3samples_0.35Induced_svgs_set",set,".csv"),header=T)
  }
  a_ = a_[,-1]
  if(kdonor ==1){a_c =a_}else{
    a_c = rbind(a_c,a_)
  }
  
}
sd_cols = colSds(data.matrix(a_c))
# Factors_df.to_csv(path.join(dir_output,"factors_nb_szMean_sample_s"+str(k+1)+"_L10_3samples_0.35Induced_svgs_set4.csv"))
#######################################
make_plot<-function(pos,
                    layer_sample_tmp___,title_,group.colors){
  # plots_l = list()
  df_tmp=data.frame( imagerow=-pos[,2],
                     imagecol= pos[,1],
                     # fill_tmp=factor_mat[,i],
                     layer=layer_sample_tmp___)
  p=ggplot(df_tmp,aes(x=imagecol,y=imagerow,fill=layer)) +
    geom_point(shape = 21, colour = "black", size = .55, stroke = NA)+
    coord_cartesian(expand=FALSE)+
    xlab("") +
    ylab("") +
    
    ggtitle(title_)+
    
    labs(fill = paste0(" "))+
    theme_set(theme_bw(base_size = 10))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size=7))+
  
    
    scale_fill_manual(values=group.colors) +
    coord_equal()
  p
}

##################
for(kdonor in c(1:3)){
  # for(k in 1){
    ksample = (kdonor-1)*4+ set 
    layers=list_layer[[ksample]]
    print(ksample)
    ##
    if(set==1){    
      a_=read.csv(paste0("/dcs04/hansen/data/ywang/ST/DLPFC/threeSamples/factors_nb_szMean_sample_s",
                         kdonor,"_L10_3samples_0.35Induced_svgs.csv"),header=T)#factors for 3-sample mNSF
    }else{
      a_=read.csv(paste0("//dcs04/hansen/data/ywang/ST/DLPFC/threeSamples/factors_nb_szMean_sample_s",
                         kdonor,"_L10_3samples_0.35Induced_svgs_set",set,".csv"),header=T)
    }
    
    a_ = a_[,-1]
    X=read.csv(paste0('//dcs04/hansen/data/ywang/ST/DLPFC/processed_Data///X_allSpots_sample',ksample,'.csv'))
    
    
    layers_pred_sampleTmp = walktrap_clustering(clusternum=ifelse(kdonor == 2, 5, 7), latent_dat=t((a_[,]))/sd_cols, knearest=70)
    which_na = which(is.na(layers))
    print(ARI(layers_pred_sampleTmp[-which_na], layers[-which_na]))
    
    
    saveRDS(layers_pred_sampleTmp, file = paste0("layers_pred_withinEachsample_sample_",ksample,"_scaleV2_3samples_010induced_set",set,".rds"))
    
  
}


#######################################

ARI_vec=c()
# dir.create("domain_mnsf_dlpfc")
# dir.create("domain_true_dlpfc")

for(kdonor in c(1:3)){
  # for(k in 1){
    ksample = (kdonor-1)*4+ set 
    layers=list_layer[[ksample]]
    
    ##
    
    X=read.csv(paste0('//dcs04/hansen/data/ywang/ST/DLPFC/processed_Data///X_allSpots_sample',ksample,'.csv'))
    
    
    layers_pred_sampleTmp = readRDS(paste0("layers_pred_withinEachsample_sample_",ksample,"_scaleV2_3samples_010induced_set",set,".rds"))
    # saveRDS(layers_pred_sampleTmp, file = paste0("layers_pred_withinEachsample_sample_",ksample,"_scaleV2_3samples_010induced.rds"))
    
    clusterlabel_refine = refine_cluster_10x(clusterlabels=layers_pred_sampleTmp,
                                             location=location_list[[ksample]],shape="hexagon")
    
    which_na = which(is.na(layers))
    
    ART_tmp = ARI(clusterlabel_refine[-which_na], layers[-which_na])
    print(ART_tmp)
    ARI_vec = c(ARI_vec, ART_tmp)
    
    saveRDS(clusterlabel_refine, file = paste0("layers_pred_withinEachsample_sample_",ksample,"_refined_scaleV2_3samples_010induced_set",set,".rds"))
    
    pdf(paste0("domain_true_dlpfc/set",set,"_donor",kdonor,".pdf"),height=2.5,width=2.5)
    
    p1= make_plot(X[,],layers[], paste0("true layer, donor ",kdonor ),
                  group.colors <- c(Layer1 = "#FC8D62", Layer2 = "#FFD92F", Layer3 ="#A6D854", Layer4 = "#66C2A5", Layer5 = "#00A9FF",
                                    Layer6="#8DA0CB",WM="#E78AC3"))
    
    print(p1 + theme(legend.position = "none"))
    dev.off()
    
    
    pdf(paste0("domain_mnsf_dlpfc/set",set,"_donor",kdonor,".pdf"),height=2.5,width=2.5)
    
    p1= make_plot(X[,], as.character(clusterlabel_refine), paste0("predicted domain by mNSF, donor ",kdonor ),
                  group.colors=c("1" = "#FC8D62", "2" = "#FFD92F", "3" ="#A6D854", "4" = "#66C2A5", "5" = "#00A9FF",
                                 "6"="#8DA0CB", "7"="#E78AC3")    )
    print(p1 + theme(legend.position = "none"))
    dev.off()
  # }
}
saveRDS(ARI_vec, file = paste0("ARI_vec_sample1_5_9_predLayers_mnsf_ClusterWithinEachSample_mnsfrefined_scaleV2_3samples_034induced_set",set,".rds"))
# set1
# 0.419376
# 0.3557843
# 0.4115804


# set2
# 0.4635874
# 0.3255512
# 0.3006927

# set3
# 0.3504535
# 0.5966603
# 0.2188306

# set 4
# 0.29984
# 0.5262868
# 0.3164823


