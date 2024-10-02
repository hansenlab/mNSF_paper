# ml conda_R/4.4
library(SpatialPCA)

### load pacakges 
library(tidyverse)
library(ggplot2)
library(Matrix)
# library(Rmisc)#
library(ggforce)#
# library(rjson)#
library(cowplot)
library(RColorBrewer)
library(grid)
require(nnet)
# require(nnet)
library(aricode)


#######
library(RColorBrewer)
library(ggplot2)
# library(readbitmap)#
myPalette_ <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
# https://ggplot2-book.org/scale-colour.html
myPalette = scale_fill_brewer(palette = "Set2")

#######
dir_Data = "/dcs04/hansen/data/ywang/revision/SpatialPCA_example_data/"
dir_dlpfc_data = paste0(dir_Data,"DLPFC/")

dir_output="/dcs04/hansen/data/ywang/ST/DLPFC/threeSamples/"
setwd(dir_output)
dir_processedData = "/dcs04/hansen/data/ywang/ST/DLPFC/processed_Data_keepAll/"
# load(paste0(dir_dlpfc_data, "LIBD_sample",ksample,".RData"))

# colnames(count_sub) = rownames(xy_coor
#######################################

#######################################
list_layer=list()
####################################### 
for(ksample in 1:12){
  
  load(paste0("//dcs04/hansen/data/ywang/archive/scRNA/Oct5_2021_Lukas_data_more_Genes/out/layer_sample_Sample",ksample,".RData"))
  
  
  list_layer[[ksample]]=layer
  
  
}

make_plot<-function(pos,
                    layer_sample_tmp___,title_,group.colors){
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

#######################################
set=1
# set=2
# set=3
# set=4
# 
out_SpatialPCA = readRDS(paste0("/dcs04/hansen/data/ywang/revision/domain/out_SpatialPCA_3samples_svgs_Set",set,"_L10.rds"))


# pdf(paste0("sample1_5_9_predLayers_mnsf_ClusterWithinEachSample_spatialPCA_svgs_set",set,".pdf"),height=2.5,width=2.5)
for(kdonor in c(1:3)){
  ksample = (kdonor-1)*4+ set
  layers=list_layer[[ksample]]
  
  ##
  a_ = t(out_SpatialPCA$SpatialPC_list[[kdonor]][1:10,])
  
  X=read.csv(paste0('//dcs04/hansen/data/ywang/ST/DLPFC/processed_Data///X_allSpots_sample',ksample,'.csv'))
  # 
  layers_pred_sampleTmp = walktrap_clustering(clusternum=ifelse(kdonor == 2, 5, 7), latent_dat=t((a_)), knearest=70)
  saveRDS(layers_pred_sampleTmp, file = paste0("layers_pred_withinEachsample_donor_",kdonor,"_spatialPCA_svgs_set",set,".rds"))

  
  which_na = which(is.na(layers))
  print(ARI(layers_pred_sampleTmp[-which_na], layers[-which_na]))
  
}



######################
dir.create("domain_spatialPCA_dlpfc")

ARI_vec=c()



pdf(paste0("sample1_5_9_predLayers_mnsf_ClusterWithinEachSample_spatialPCA_refined_svgs_set",set,".pdf"),height=2.5,width=2.5)
for(kdonor in c(1:3)){
  ksample = (kdonor-1)*4+ set
  layers=list_layer[[ksample]]
  
  ##
  a_ = t(out_SpatialPCA$SpatialPC_list[[kdonor]][1:10,])

  X=read.csv(paste0('//dcs04/hansen/data/ywang/ST/DLPFC/processed_Data///X_allSpots_sample',ksample,'.csv'))

  layers_pred_sampleTmp = readRDS(file = paste0("layers_pred_withinEachsample_donor_",kdonor,"_spatialPCA_svgs_set",set,".rds"))
  layers_pred_sampleTmp_Refined = refine_cluster_10x(clusterlabels=layers_pred_sampleTmp,location=X,shape="hexagon")
  
  pdf(paste0("domain_spatialPCA_dlpfc/set",set,"_donor",kdonor,".pdf"),height=2.5,width=2.5)
  

  p1= make_plot(X[,], as.character(layers_pred_sampleTmp_Refined), paste0("predicted domain by SpatialPCA, donor ",kdonor ),
                group.colors=c("1" = "#FC8D62", "2" = "#FFD92F", "3" ="#A6D854", "4" = "#66C2A5", "5" = "#00A9FF",
                               "6"="#8DA0CB", "7"="#E78AC3")    )
  print(p1 + theme(legend.position = "none"))
  dev.off()
  
  
  
  saveRDS(layers_pred_sampleTmp_Refined, file = paste0("layers_pred_withinEachsample_donor_",kdonor,"_spatialPCA_refined_svgs_set",set,".rds"))
  which_na = which(is.na(layers))
  ART_tmp = ARI(layers_pred_sampleTmp_Refined[-which_na], layers[-which_na])
  print(ART_tmp)
  ARI_vec = c(ARI_vec, ART_tmp)
  saveRDS(layers_pred_sampleTmp_Refined, file = paste0("layers_pred_withinEachsample_donor_",kdonor,"_spatialPCA_refined_svgs_",set,".rds"))
  
  
}
saveRDS(ARI_vec,file = paste0("ARI_vec_sample1_5_9_predLayers_mnsf_ClusterWithinEachSample_spatialPCA_refined_svgs_set",set,"_L10.rds"))
# set1

# set2：
# [1] 0.2962618
# [1] 0.4175159
# [1] 0.3272038
# set3:
# 0.3815678
# 0.3405254
# 0.3140999
# set4:
# 0.3269615
# 0.3950879
# 0.3162719
